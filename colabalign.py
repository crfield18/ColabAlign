'''Top-level driver for structure alignment, clustering, and consensus selection.'''

from os import cpu_count
import warnings
from pathlib import Path
from time import time

from cli import script_args
from structure_prep import discover_input_models, prepare_models
from usalign_runner import run_pairwise_alignments, check_tm_matrix_completeness
from dendrogram import grow_tree as _grow_tree, tree_clustering as _tree_clustering
from consensus import ClusterConsensus


class ColabAlign:
    '''Coordinates structure prep, pairwise alignment, clustering, and consensus selection.'''

    def __init__(self, args) -> None:
        '''Initialize output directories and runtime parameters from parsed CLI args.'''
        self.model_list = discover_input_models(args.input)
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

        self.usalign_path = args.usalign
        self.beem_path = args.beem
        self.mode = args.mode

        self.tmmatrix = None
        self.library = None
        self.transforms = None
        self.clusters_by_threshold = None

    def prepare_structures(self):
        '''Prepare input models (e.g. chain extraction, format conversion) for alignment.'''
        self.model_list = prepare_models(
            self.model_list, self.models_path, self.beem_path, self.mode, self.cores
            )

    def pairwise_alignment(self):
        '''Run all-vs-all USalign alignments and write the resulting TM-score matrix.'''
        cache_dir = self.results_path.joinpath('usalign_cache')
        self.tmmatrix, self.library, self.transforms, failed_pairs = run_pairwise_alignments(
            self.model_list, self.models_path, self.usalign_path, self.cores, cache_dir
        )
        self.tmmatrix.to_csv(self.results_path.joinpath('us-align_score_matrix.csv'), index=True)

        missing = check_tm_matrix_completeness(self.tmmatrix)
        if missing:
            total_pairs = len(self.model_list) * (len(self.model_list) - 1) // 2
            print(
                f'\nWARNING: {len(missing)} of {total_pairs} pairs are missing from the '
                f'TM-score matrix (see failed-pair list above). The dendrogram and '
                f'clustering below will proceed, but any structures only connected '
                f'through a missing pair may cluster unexpectedly.'
            )

    def grow_tree(self):
        '''Build and save the dendrogram from the TM-score matrix.'''
        _grow_tree(self.tmmatrix, self.model_list, self.results_path)

    def tree_clustering(self):
        '''Cluster the dendrogram at each configured threshold.'''
        self.clusters_by_threshold = _tree_clustering(
            self.thresholds, self.results_path, self.tmmatrix,
            self.transforms, self.model_list, self.models_path
        )

    def find_representatives(self):
        '''Select and write representative structures for each cluster at each threshold.'''
        reps_parser = ClusterConsensus(
            library=self.library, tmmatrix=self.tmmatrix, threads=self.cores
        )
        reps_parser.get_representatives_all_thresholds(
            self.results_path, self.clusters_by_threshold
            )


def main():
    '''Run the full ColabAlign pipeline end-to-end and report timing for each stage.'''
    main_time_start = time()
    local_instance = ColabAlign(script_args())

    prep_time_start = time()
    print('Preparing input structures.')
    local_instance.prepare_structures()
    print(f'Structure preparation elapsed time:\t{time() - prep_time_start:.2f} s')

    align_time_start = time()
    print('Starting pairwise alignment.')
    local_instance.pairwise_alignment()
    print(f'Pairwise alignment elapsed time:\t{time() - align_time_start:.2f} s')

    tree_time_start = time()
    print('Starting dendrogram calculation.')
    local_instance.grow_tree()
    print(f'Dendrogram calculation elapsed time:\t{time() - tree_time_start:.2f} s')

    tree_cluster_time_start = time()
    print('Starting dendrogram clustering.')
    local_instance.tree_clustering()
    print(f'Dendrogram clustering elapsed time:\t{time() - tree_cluster_time_start:.2f} s')

    reps_time_start = time()
    print('Finding cluster representatives.')
    local_instance.find_representatives()
    print(f'Representative calculation elapsed time:\t{time() - reps_time_start:.2f} s')

    print(f'Total elapsed time:\t{time() - main_time_start:.2f} s')


if __name__ == '__main__':
    main()
