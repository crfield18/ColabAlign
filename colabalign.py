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
    def _reverse_transformation_matrix(self, transform_mx_forward):
        assert isinstance(transform_mx_forward, np.ndarray)
        assert transform_mx_forward.shape == (3, 4)

        translate_v_forward = transform_mx_forward[:, 0:1]
        rotation_mx_forward = transform_mx_forward[:, 1:4]

        rotation_mx_reverse = rotation_mx_forward.T

        # This step is matrix multiplication NOT a dot product (np.dot handles both)
        translate_v_reverse = - np.dot(rotation_mx_reverse, translate_v_forward)
        transform_mx_reverse = np.hstack((translate_v_reverse, rotation_mx_reverse))

        return transform_mx_reverse

    def _run_usalign(self, model1, model2):
        # Align model against models listed in corresponding model list file
        cmd = [
            self.usalign_path,
            self.models_path.joinpath(model1.name).as_posix(),
            self.models_path.joinpath(model2.name).as_posix(),
            '-outfmt', '2', '-m', '-'
        ]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            stdout, stderr = process.communicate()
        return model1, model2, stdout, stderr

    def _usalign_process(self, job):
        results = []
        for pair in job:
            results.append(self._run_usalign(pair[0], pair[1]))
        return results

    def _parse_usalign_stdout(self, stdout:str):
        decoded_results = stdout.decode('UTF8').split('\n')
        alignment = decoded_results[1].split('\t')
        tm_1 = alignment[2]
        tm_2 = alignment[3]

        transform_mx_forward = np.array([[d for d in line.split(' ') if d][1:]
                                         for line in decoded_results[4:7]]).astype(np.float64)
        transform_mx_reverse = self._reverse_transformation_matrix(
            transform_mx_forward.astype(np.float64)
        )
        return tm_1, tm_2, transform_mx_forward, transform_mx_reverse

    def _results_map_to_df(self, results_dict: dict) -> pd.DataFrame:
        def generate_tuples():
            for model1, model2_dict in results_dict.items():
                for model2, value in model2_dict.items():
                    yield (model1, model2, value['tm_forward'],
                           value['transform_mx_forward'],
                           value['tm_reverse'],
                           value['transform_mx_reverse'])

        return pd.DataFrame(generate_tuples(),
                            columns=['model1', 'model2', 'tm_forward',
                                     'transform_mx_forward', 'tm_reverse',
                                     'transform_mx_reverse'])

    def _find_matrix_in_array(self, row, model):
        if row['model1'] == model:
            return row['transform_mx_forward']
        if row['model2'] == model:
            return row['transform_mx_reverse']
        return None

    # Make structural dendrogram from US-align results
    def pairwise_alignment(self):
        def _chunks(lst, n):
            '''Yield n chunks from the list.'''
            avg = len(lst) / n
            last = 0
            while last < len(lst):
                yield lst[int(last):int(last + avg)]
                last += avg

        def _run_beem(self, cif_file):
            cmd = [
                self.beem_path,
                self.models_path.joinpath(cif_file).as_posix(),
                '-outfmt=2',
                f'-p={cif_file.stem}'
            ]
            with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
                stdout, stderr = process.communicate()
            return stdout, stderr

        print('Copying pdb files and converting .cif files to .pdb.')

        aa_codes = {
            'CYS', 'ASP', 'SER', 'GLN', 'LYS', 'ILE', 'PRO', 'THR',
            'PHE', 'ASN', 'GLY', 'HIS', 'LEU', 'ARG', 'TRP', 'ALA',
            'VAL', 'GLU', 'TYR', 'MET'
        }

        original_stems = {m.stem for m in self.model_list}
        for model in self.model_list:
            # Convert any .cif files to .pdb format for MUSTANG compatibility.
            # BeEM also splits multi-model .cif files into individual .pdb files with chain IDs appended to the end of the filename.
            # We achieve the same for native pdb files using biopython.
            if model.suffix == '.cif':
                beem_output, _ = _run_beem(self, model)
                beem_models = (item for item in beem_output.decode(errors='ignore').split('\n') if item != '')
                for m in beem_models:
                    shutil.move(src=m, dst=self.models_path.joinpath(m))

            else:
                parser = PDBParser(QUIET=True)
                io = PDBIO()
                structure = parser.get_structure(model.stem, model)

                for chain in structure.get_chains():
                    chain_id = chain.get_id().upper()
                    out_name = f'{structure.get_id()}{chain_id}.pdb'
                    io.set_structure(chain)
                    io.save((self.models_path.joinpath(out_name)).as_posix())

        # Delete any files that are empty or do not contain protein residues as
        # US-align does not handle them properly.
        for structure_file in self.models_path.glob('*.pdb'):
            if not structure_file.stat().st_size > 0 or not structure_file.is_file():
                print(f'Removing empty file: {structure_file}')
                Path.unlink(structure_file)

            with open(structure_file, 'r', encoding='utf-8') as f:
                file_contains_protein = False
                for line in f:
                    if line.startswith('ATOM') and line[17:20].strip() in aa_codes:
                        file_contains_protein = True
                        break

                if not file_contains_protein:
                    print(f'Removing non-protein file: {structure_file}')
                    Path.unlink(structure_file)
                else:
                    continue

        self.model_list = sorted([f for f in self.models_path.glob('*.pdb')])

        # Only take the first chain (alphabetically, most of the time it's chain A)
        # if mode is set to first
        if self.mode == 'first':
            seen_inputs = set()
            first_chains = []
            for model in self.model_list:
                input_stem = model.stem[:-1]
                if input_stem in original_stems and input_stem not in seen_inputs:
                    seen_inputs.add(input_stem)
                    first_chains.append(model)
                else:
                    model.unlink()
            self.model_list = first_chains

        # Calculate all possible combinations of models, then create discrete lists of comparisons
        # on a per-model basis. This enables us to run multiple, concurrent US-align instances and
        # avoiding doing any redundant calculations
        all_combos = combinations(self.model_list, 2)

        # Run multiple US-align jobs in parallel
        # This greatly speeds up the workflow
        print(f'Running US-align jobs for {len(self.model_list)} chains.')
        process_results = []  # Collect dictionaries from each future

        with ProcessPoolExecutor(max_workers=self.cores) as executor:
            # Build jobs (chunks of pairwise combos)
            combos_list = list(all_combos)

            # Splitting the jobs into a lot of chunks to make the progess bar have any use at all
            if len(combos_list) <= int(self.cores*500):
                chunk_count = int(self.cores)
            else:
                chunk_count = int(self.cores*500)

            jobs = list(_chunks(combos_list, chunk_count))

            total_bar = tqdm_auto(
                total=len(combos_list),
                desc='Total alignments',
                unit=' alignment',
                dynamic_ncols=True,
                leave=True,
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
