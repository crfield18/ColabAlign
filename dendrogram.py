'''Dendrogram construction, tree clustering, and per-cluster structural superposition.'''

from math import isnan
from pathlib import Path
from shutil import copy
import subprocess

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.PDB import PDBParser, PDBIO

from usalign_runner import get_transform_matrix


class StructureAligner:  # pylint: disable=too-few-public-methods
    '''Applies a USalign rotation/translation transform to a structure's coordinates.'''

    def __init__(self, input_file: Path, output_file: Path, transform_matrix: np.ndarray) -> None:
        '''Validate inputs and store the transform to apply.'''
        self.input_file = input_file
        self.output_file = output_file
        self.transform_matrix = transform_matrix.astype(np.float64)
        assert transform_matrix.shape == (3, 4)

        if input_file.suffix not in ('.pdb', '.cif'):
            raise ValueError(f'Invalid file extension: {input_file.suffix}. Expected .pdb or .cif')

    @staticmethod
    def _apply_transform(atom, rotation_matrix, translation_vector):
        '''Apply the rotation/translation transform to a single atom's coordinates.'''
        new_coord = np.dot(rotation_matrix, atom.get_coord()) + translation_vector.flatten()
        atom.set_coord(new_coord)

    def _transform_atom(self, atom, rotation_matrix, translation_vector):
        '''Apply the transform to an atom, handling disordered (altloc) atoms.'''
        if atom.is_disordered():
            for altloc in atom.disordered_get_id_list():
                child = atom.disordered_get(altloc)
                self._apply_transform(child, rotation_matrix, translation_vector)
        else:
            self._apply_transform(atom, rotation_matrix, translation_vector)

    def transform_coords(self) -> Path:
        '''Parse the input structure, apply the transform to every atom, and save it.'''
        parser = PDBParser(QUIET=True)
        io = PDBIO()
        structure = parser.get_structure('structure', self.input_file)

        rotation_matrix = self.transform_matrix[:, 1:4].reshape((3, 3))
        translation_vector = self.transform_matrix[:, 0].reshape((1, 3))

        for atom in structure.get_atoms():
            self._transform_atom(atom, rotation_matrix, translation_vector)

        io.set_structure(structure)
        io.save(str(self.output_file))
        return self.output_file


def grow_tree(tmmatrix: pd.DataFrame, model_list, results_path: Path):
    '''Build a UPGMA dendrogram from the TM-score matrix and write it as Newick.'''
    print('Generating structural dendrogram.')

    distances_df = 1.0000 - tmmatrix
    distances_df = distances_df.round(4)
    distances_df = distances_df.clip(lower=0.0000, upper=1.0000)

    lower_tri_df = distances_df.where(np.tril(np.ones(distances_df.shape)).astype(bool))
    lower_tri_lists = [
        [value for value in row if not isnan(value)]
        for row in lower_tri_df.values.tolist()
    ]

    tm_matrix = DistanceMatrix(
        names=[model.stem for model in model_list],
        matrix=lower_tri_lists
    )
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(tm_matrix)

    Phylo.write(tree, results_path.joinpath('colabalign.tree'), 'newick')


def _run_treecluster(thresholds, results_path: Path):
    '''Run TreeCluster.py for each threshold, in parallel, and wait for completion.'''
    completed_processes = []
    for thresh in thresholds:
        cmd = [
            'TreeCluster.py',
            '-i', results_path.joinpath('colabalign.tree').as_posix(),
            '-o', results_path.joinpath(f'clusters/{thresh:.2f}.tsv').as_posix(),
            '-t', f'{thresh:.2f}'
        ]
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            process.communicate()
            completed_processes.append(process)
    for process in completed_processes:
        process.wait()


def _parse_cluster_file(results_path: Path, threshold_str: str):
    '''Parse a TreeCluster.py output TSV into dict[cluster_num -> list of member names].'''
    clusters = {}
    cluster_file_path = results_path.joinpath(f'clusters/{threshold_str}.tsv')
    with open(cluster_file_path, 'r', encoding='UTF8') as cluster_file:
        for line in cluster_file:
            if line.startswith('SequenceName'):
                continue
            name, cluster_num = line.strip('\n').split('\t')
            clusters.setdefault(int(cluster_num), []).append(name.strip())
    return clusters


def _copy_unclustered_members(member_names, model_list, models_path, cluster_output_dir):
    '''Copy singleton/unclustered (-1) members over without any superposition.'''
    for name in member_names:
        original_file = next((f for f in model_list if f.stem == name), None)
        if original_file:
            copy(models_path.joinpath(original_file.name),
                 cluster_output_dir.joinpath(original_file.name))


def _select_medoid(tmmatrix: pd.DataFrame, member_names):
    '''Pick the member with the highest mean TM-score to the rest of the cluster.'''
    sub_matrix = tmmatrix.loc[list(member_names), list(member_names)]
    return sub_matrix.mean(axis=0).idxmax()


def _superpose_member(
    name, reference_model_name, transforms, model_list, models_path, cluster_output_dir
    ):
    '''Superpose a single cluster member onto the reference (medoid) structure.'''
    transform_matrix = get_transform_matrix(transforms, name, reference_model_name)
    if transform_matrix is None:
        print(
            f'Warning: no transform found for {name} -> {reference_model_name}, skipping superposition'
            )
        return

    model_file = next((f for f in model_list if f.stem == name), None)
    if model_file is None:
        return

    aligner = StructureAligner(
        input_file=models_path.joinpath(model_file.name),
        output_file=cluster_output_dir.joinpath(model_file.name),
        transform_matrix=transform_matrix
    )
    aligner.transform_coords()


def _process_cluster(cluster_num, member_names, ctx, cluster_output_dir):
    '''Copy the medoid and superpose the rest of a single cluster's members onto it.'''
    tmmatrix, transforms, model_list, models_path = ctx

    if cluster_num == -1:
        _copy_unclustered_members(member_names, model_list, models_path, cluster_output_dir)
        return

    # Medoid selection: highest MEAN TM-score to the rest of the cluster,
    # not the single highest pairwise TM-score to any one member (which
    # can pick an atypical structure).
    reference_model_name = _select_medoid(tmmatrix, member_names)

    ref_file = next((f for f in model_list if f.stem == reference_model_name), None)
    if ref_file is None:
        return

    copy(models_path.joinpath(ref_file.name), cluster_output_dir.joinpath(ref_file.name))

    for name in member_names:
        if name == reference_model_name:
            continue
        _superpose_member(
            name, reference_model_name, transforms, model_list, models_path, cluster_output_dir
            )


# pylint: disable=too-many-arguments,too-many-positional-arguments
def tree_clustering(thresholds, results_path: Path, tmmatrix: pd.DataFrame,
                     transforms: dict, model_list, models_path: Path):
    '''
    Runs TreeCluster.py per threshold, superposes cluster members onto the
    cluster's medoid (highest MEAN TM-score to the rest of the cluster,
    not just the single highest pairwise TM-score), and returns
    dict[threshold_str -> dict[cluster_num -> list of member stems]].
    '''
    print('Calculating clusters from structural tree.')

    _run_treecluster(thresholds, results_path)

    ctx = (tmmatrix, transforms, model_list, models_path)
    all_clusters_by_threshold = {}

    for thresh in thresholds:
        threshold_str = f'{thresh:.2f}'
        clusters = _parse_cluster_file(results_path, threshold_str)

        for cluster_num, member_names in clusters.items():
            cluster_output_dir = results_path.joinpath(f'clusters/{threshold_str}/{cluster_num}')
            cluster_output_dir.mkdir(exist_ok=True, parents=True)
            _process_cluster(cluster_num, member_names, ctx, cluster_output_dir)

        all_clusters_by_threshold[threshold_str] = clusters

    return all_clusters_by_threshold
