'''Top-level driver for structure alignment, clustering, and consensus selection.'''

from os import cpu_count
import shutil
import warnings
from types import SimpleNamespace
from pathlib import Path
from time import time

from .cli import script_args
from .structure_prep import discover_input_models, prepare_models
from .usalign_runner import run_pairwise_alignments, check_tm_matrix_completeness
from .dendrogram import grow_tree as _grow_tree, tree_clustering as _tree_clustering
from .consensus import ClusterConsensus


def _require_executable(path: Path, tool: str, flag):
    '''Resolve an executable on PATH (or as given), failing early with a clear message.'''
    resolved = shutil.which(str(path))
    if resolved is None:
        hint = f' or point to it with {flag}' if flag else ''
        raise SystemExit(
            f'Could not find {tool} ({path}). Install it (conda install -c bioconda '
            f'colabalign pulls in all external tools){hint}.'
        )
    return Path(resolved)


class ColabAlign:
    '''Coordinates structure prep, pairwise alignment, clustering, and consensus selection.'''

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def __init__(self, inputs, output, cores: int = 1, threshold=(0.20,),
                 usalign='USalign', beem='BeEM', mode: str = 'all') -> None:
        '''
        Set up output directories and runtime parameters.

        inputs: iterable of .pdb/.cif files, directories, or .zip archives.
        output: output directory (created if missing).
        threshold: one clustering threshold or an iterable of them.
        mode: 'all' chains or only the 'first' chain of each input.
        '''
        if mode not in ('all', 'first'):
            raise ValueError(f"mode must be 'all' or 'first', not {mode!r}")
        if isinstance(inputs, (str, Path)):
            inputs = [inputs]
        if isinstance(threshold, (int, float)):
            threshold = [threshold]
        args = SimpleNamespace(
            input=[Path(p) for p in inputs], output=Path(output), cores=int(cores),
            threshold=[float(t) for t in threshold], usalign=Path(usalign),
            beem=Path(beem), mode=mode,
        )

        self.thresholds = args.threshold

        print('Setting up file structure.')
        self.home_path = args.output
        self.models_path = self.home_path.joinpath('models')
        self.results_path = self.home_path.joinpath('results')

        for output_subdir in (
            self.home_path, self.models_path, self.results_path,
            self.results_path.joinpath('clusters'),
        ):
            Path.mkdir(output_subdir, exist_ok=True, parents=True)

        # Temporary home for the contents of any uploaded .zip archives;
        # removed once structure preparation has finished with them.
        self.extract_path = self.home_path.joinpath('.extracted_inputs')
        self.model_list = discover_input_models(args.input, self.extract_path)
        if not self.model_list:
            raise SystemExit('No .pdb, .cif or .zip input files found.')

        self.cores = args.cores
        if self.cores <= 0:
            warnings.warn(
                'User-defined core count cannot be lower than 1. '
                'Core count set to default (1).'
            )
            self.cores = 1
        if self.cores > cpu_count():
            warnings.warn(
                f'User-defined core count ({self.cores}) exceeds available cores '
                f'({cpu_count()}). Using maximum available cores instead.'
            )
            self.cores = cpu_count()

        self.usalign_path = _require_executable(args.usalign, 'USalign', '-u/--usalign')
        if any(m.suffix.lower() == '.cif' for m in self.model_list):
            self.beem_path = _require_executable(args.beem, 'BeEM', '-b/--beem')
        else:
            self.beem_path = args.beem  # only needed for .cif inputs
        _require_executable(Path('TreeCluster.py'), 'TreeCluster.py', None)
        self.mode = args.mode

        self.tmmatrix = None
        self.pair_store = None
        self.clusters_by_threshold = None

    def prepare_structures(self):
        '''Prepare input models (e.g. chain extraction, format conversion) for alignment.'''
        self.model_list = prepare_models(
            self.model_list, self.models_path, self.beem_path, self.mode, self.cores)
        shutil.rmtree(self.extract_path, ignore_errors=True)

    def pairwise_alignment(self):
        '''Run all-vs-all USalign alignments and write the resulting TM-score matrix.'''
        cache_dir = self.results_path.joinpath('usalign_cache')
        self.tmmatrix, self.pair_store, _ = run_pairwise_alignments(
            self.model_list, self.usalign_path, self.cores, cache_dir
        )
        self.tmmatrix.to_csv(self.results_path.joinpath('us-align_score_matrix.csv'), index=True)

        missing = check_tm_matrix_completeness(self.tmmatrix)
        if missing:
            total_pairs = len(self.model_list) * (len(self.model_list) - 1) // 2
            print(
                f'\nWARNING: {len(missing)} of {total_pairs} pairs are missing from the '
                f'TM-score matrix (see failed-pair list above). Missing pairs are '
                f'treated as unrelated (TM-score 0) for the dendrogram and clustering, '
                f'so any structures only connected through a missing pair may cluster '
                f'unexpectedly.'
            )
            # NaNs break the lower-triangle DistanceMatrix and the per-cluster
            # guide trees. The CSV written above keeps them, for transparency.
            self.tmmatrix = self.tmmatrix.fillna(0.0)

    def grow_tree(self):
        '''Build and save the dendrogram from the TM-score matrix.'''
        _grow_tree(self.tmmatrix, self.model_list, self.results_path)

    def tree_clustering(self):
        '''Cluster the dendrogram at each configured threshold.'''
        self.clusters_by_threshold = _tree_clustering(
            self.thresholds, self.results_path, self.tmmatrix,
            self.pair_store, self.model_list, self.models_path
        )

    def find_representatives(self):
        '''Select and write representative structures for each cluster at each threshold.'''
        reps_parser = ClusterConsensus(
            pair_store=self.pair_store, tmmatrix=self.tmmatrix, threads=self.cores
        )
        reps_parser.get_representatives_all_thresholds(
            self.results_path, self.clusters_by_threshold
            )


    @classmethod
    def from_args(cls, args):
        '''Build an instance from an argparse Namespace produced by cli.script_args().'''
        return cls(inputs=args.input, output=args.output, cores=args.cores,
                   threshold=args.threshold, usalign=args.usalign, beem=args.beem,
                   mode=args.mode)

    def run(self):
        '''Run the full pipeline end-to-end, reporting the elapsed time of each stage.'''
        main_time_start = time()
        stages = (
            ('Preparing input structures.', self.prepare_structures,
             'Structure preparation'),
            ('Starting pairwise alignment.', self.pairwise_alignment,
             'Pairwise alignment'),
            ('Starting dendrogram calculation.', self.grow_tree,
             'Dendrogram calculation'),
            ('Starting dendrogram clustering.', self.tree_clustering,
             'Dendrogram clustering'),
            ('Finding cluster representatives.', self.find_representatives,
             'Representative calculation'),
        )
        for start_message, stage, label in stages:
            stage_time_start = time()
            print(start_message)
            stage()
            print(f'{label} elapsed time:\t{time() - stage_time_start:.2f} s')
        print(f'Total elapsed time:\t{time() - main_time_start:.2f} s')


def main(argv=None):
    '''Entry point for the colabalign command.'''
    ColabAlign.from_args(script_args(argv)).run()


if __name__ == '__main__':
    main()
