'''Top-level driver for structure alignment, clustering, and consensus selection.'''

from os import cpu_count
import shutil
import warnings
from types import SimpleNamespace
from pathlib import Path
from time import time

from .cli import script_args
from .structure_prep import (
    collect_input_models, expand_duplicate_rows, prepare_models_with_chains, write_duplicates_csv,
)
from .usalign_runner import run_pairwise_alignments, check_tm_matrix_completeness
from .dendrogram import grow_tree as _grow_tree, tree_clustering as _tree_clustering
from .consensus import ClusterConsensus
from .panel import (
    DEFAULT_LOW_COVERAGE_FRACTION, MAX_LOW_COVERAGE_FRACTION, build_panel, choose_cuts,
    cut_options, format_threshold, parse_panel_size, write_cut_options, write_panel_summary,
)

DEFAULT_THRESHOLD = 0.20


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


def _normalise_panel_sizes(panel_size):
    '''None, 50, "40-60", (40, 60) or a list of these -> list of (min, max) tuples.'''
    if panel_size is None:
        return []
    if isinstance(panel_size, (int, str)) or (
            isinstance(panel_size, tuple) and all(isinstance(v, int) for v in panel_size)):
        panel_size = [panel_size]
    return [parse_panel_size('-'.join(map(str, size)) if isinstance(size, tuple) else size)
            for size in panel_size]


class ColabAlign:
    '''Coordinates structure prep, pairwise alignment, clustering, and consensus selection.'''

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def __init__(self, inputs, output, cores: int = 1, threshold=None,
                 usalign='USalign', beem='BeEM', mode: str = 'all', panel_size=None,
                 low_coverage_fraction: float = DEFAULT_LOW_COVERAGE_FRACTION) -> None:
        '''
        Set up output directories and runtime parameters.

        inputs: iterable of .pdb/.cif files, directories, or .zip/.tar(.gz) archives.
        output: output directory (created if missing).
        threshold: one clustering threshold or an iterable of them. Defaults to
            0.20 when no panel_size is given.
        panel_size: desired panel size(s); each an int, a (min, max) tuple or a
            "MIN-MAX" string. A threshold is chosen for each.
        low_coverage_fraction: fraction of each panel given to the structures
            least covered by the cluster representatives.
        mode: 'all' chains or only the 'first' chain of each input.
        '''
        if mode not in ('all', 'first'):
            raise ValueError(f"mode must be 'all' or 'first', not {mode!r}")
        if not 0 <= low_coverage_fraction <= MAX_LOW_COVERAGE_FRACTION:
            raise ValueError(
                f'low_coverage_fraction must be between 0 and {MAX_LOW_COVERAGE_FRACTION}')
        if isinstance(inputs, (str, Path)):
            inputs = [inputs]
        if isinstance(threshold, (int, float)):
            threshold = [threshold]
        panel_size = _normalise_panel_sizes(panel_size)
        if threshold is None:
            threshold = [] if panel_size else [DEFAULT_THRESHOLD]
        args = SimpleNamespace(
            input=[Path(p) for p in inputs], output=Path(output), cores=int(cores),
            threshold=[float(t) for t in threshold], usalign=Path(usalign),
            beem=Path(beem), mode=mode,
        )

        self.thresholds = args.threshold
        self.panel_sizes = panel_size
        self.low_coverage_fraction = float(low_coverage_fraction)
        self.requested_by_threshold = {}
        self.expected_clusters = {}

        print('Setting up file structure.')
        self.home_path = args.output
        self.models_path = self.home_path.joinpath('models')
        self.results_path = self.home_path.joinpath('results')

        for output_subdir in (
            self.home_path, self.models_path, self.results_path,
            self.results_path.joinpath('clusters'),
        ):
            Path.mkdir(output_subdir, exist_ok=True, parents=True)

        # Temporary home for the contents of any uploaded .zip/.tar archives;
        # removed once structure preparation has finished with them.
        self.extract_path = self.home_path.joinpath('.extracted_inputs')
        self.model_list, self.duplicate_rows = collect_input_models(args.input, self.extract_path)
        if not self.model_list:
            raise SystemExit('No .pdb, .cif, .zip or .tar(.gz) input files found.')

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
        self.representatives_by_threshold = None

    def prepare_structures(self):
        '''Prepare input models (e.g. chain extraction, format conversion) for alignment.'''
        n_renamed = sum(r['name_in_colabalign'] != r['input_file_name']
                        for r in self.duplicate_rows)
        self.model_list, chains_by_input = prepare_models_with_chains(
            self.model_list, self.models_path, self.beem_path, self.mode, self.cores)
        # Written after chain splitting so names match the matrix/tree labels.
        duplicates_csv = self.results_path.joinpath('duplicates.csv')
        write_duplicates_csv(expand_duplicate_rows(self.duplicate_rows, chains_by_input),
                             duplicates_csv)
        if self.duplicate_rows:
            print(f'{n_renamed} input file(s) shared a name with another input and were '
                  f'given a _dupN name in ColabAlign. See {duplicates_csv}.')
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

    def choose_thresholds(self):
        '''Find the threshold(s) giving each requested panel size, and list every
        achievable panel size in results/panel_sizes.csv.'''
        options = cut_options(self.results_path.joinpath('colabalign.tree'),
                              self.low_coverage_fraction)
        chosen = choose_cuts(options, self.panel_sizes) if self.panel_sizes else []
        write_cut_options(options, self.results_path.joinpath('panel_sizes.csv'))
        for option in chosen:
            threshold = option.threshold
            key = format_threshold(threshold)
            self.requested_by_threshold.setdefault(key, []).extend(option.requested)
            self.expected_clusters[key] = option.n_clusters
            if all(format_threshold(t) != key for t in self.thresholds):
                self.thresholds.append(threshold)
            note = '' if all(
                parse_panel_size(r)[0] <= option.panel_size <= parse_panel_size(r)[1]
                for r in option.requested
            ) else ' (closest achievable)'
            print(f'Panel size {", ".join(option.requested)}: threshold {key} gives '
                  f'{option.n_clusters} clusters + {option.n_low_coverage} low-coverage '
                  f'= {option.panel_size} structures{note}.')
        if self.panel_sizes:
            print(f'Every achievable panel size is listed in '
                  f'{self.results_path.joinpath("panel_sizes.csv")}.')

    def tree_clustering(self):
        '''Cluster the dendrogram at each configured threshold.'''
        self.clusters_by_threshold = _tree_clustering(
            self.thresholds, self.results_path, self.tmmatrix,
            self.pair_store, self.model_list, self.models_path
        )
        for key, expected in self.expected_clusters.items():
            clusters = self.clusters_by_threshold[key]
            found = len(clusters.get(-1, [])) + sum(1 for c in clusters if c != -1)
            if found != expected:
                print(f'WARNING: threshold {key} gave {found} clusters, not the '
                      f'{expected} predicted from the tree. The panel uses the actual clusters.')

    def find_representatives(self):
        '''Select and write representative structures for each cluster at each threshold.'''
        reps_parser = ClusterConsensus(
            pair_store=self.pair_store, tmmatrix=self.tmmatrix, threads=self.cores
        )
        self.representatives_by_threshold = reps_parser.get_representatives_all_thresholds(
            self.results_path, self.clusters_by_threshold
            )

    def assemble_panels(self):
        '''Write each threshold's panel (representatives, singletons and low-coverage
        additions) and a summary of how well each panel covers every structure.'''
        models_by_name = {m.stem: self.models_path.joinpath(m.name) for m in self.model_list}
        rows = []
        for threshold_str, clusters in self.clusters_by_threshold.items():
            rows.append(build_panel(
                threshold_str, clusters, self.representatives_by_threshold.get(threshold_str, {}),
                self.tmmatrix, self.low_coverage_fraction, models_by_name,
                self.results_path.joinpath(f'clusters/{threshold_str}'),
                self.requested_by_threshold.get(threshold_str),
            ))
        rows.sort(key=lambda row: row['panel_size'])
        write_panel_summary(rows, self.results_path.joinpath('panels.csv'))
        for row in rows:
            requested = (f' (requested {row["requested_panel_size"]})'
                         if row['requested_panel_size'] else '')
            print(f'Threshold {row["threshold"]}: {row["panel_size"]} structures{requested}; '
                  f'mean best TM-score to the panel {row["mean_best_tm"]:.3f}, '
                  f'lowest {row["min_best_tm"]:.3f}. See {row["panel_csv"]}.')


    @classmethod
    def from_args(cls, args):
        '''Build an instance from an argparse Namespace produced by cli.script_args().'''
        return cls(inputs=args.input, output=args.output, cores=args.cores,
                   threshold=args.threshold, usalign=args.usalign, beem=args.beem,
                   mode=args.mode, panel_size=args.panel_size,
                   low_coverage_fraction=args.low_coverage_fraction)

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
            ('Finding achievable panel sizes.', self.choose_thresholds,
             'Panel-size sweep'),
            ('Starting dendrogram clustering.', self.tree_clustering,
             'Dendrogram clustering'),
            ('Finding cluster representatives.', self.find_representatives,
             'Representative calculation'),
            ('Assembling panels.', self.assemble_panels, 'Panel assembly'),
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
