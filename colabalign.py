from os import cpu_count
from pathlib import Path
from argparse import ArgumentParser
from math import isnan
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
import subprocess
import warnings
from shutil import copy
from collections import defaultdict

import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

def script_args():
  parser = ArgumentParser(description='Script description goes here.')

  # Required arguments
  parser.add_argument(
    '-i', '--input', 
    type=Path,
    required=True,
    nargs='+',
    help='Paths to input files.'
  )
  parser.add_argument(
    '-o', '--output', 
    type=Path,
    required=True,
    help='Path to output directory.'
  )

  # Optional arguments
  parser.add_argument(
    '-c', '--cores', 
    type=int,
    required=False,
    default=1,
    help='Number of CPU cores to use (Default = 1).'
  )

  parser.add_argument(
    '-t', '--threshold', 
    type=float,
    required=False,
    nargs='+',
    default=[0.25],
    help='Set thresholds for clustering (Default = 0.25). Multiple values can be provided.'
  )

  parser.add_argument(
    '-p', '--path', 
    type=Path,
    required=False,
    default='USalign',
    help='Path to USalign executable. Not required if using conda install.'
  )

  return parser.parse_args()

class StructureAligner():
  def __init__(self, input_pdb:Path, output_pdb:Path, transform_matrix:np.ndarray) -> None:
    self.input_pdb = input_pdb
    self.output_pdb = output_pdb
    self.transform_matrix = transform_matrix.astype(np.float64)
    assert transform_matrix.shape == (3, 4)

    if input_pdb.suffix not in ('.pdb', '.cif'):
      raise ValueError(f'Invalid file extension: {input_pdb.suffix}. Expected .pdb or .cif')

    self.file_type = input_pdb.suffix

  def transform_coords(self) -> Path:
    with (open(self.input_pdb, 'r', encoding='UTF8') as input_file,
      open(self.output_pdb, 'w', encoding='UTF8') as output_file):
      for line in input_file:
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
          model_parser = UpdateCoords(
            atom_line=line,
            file_format=self.file_type,
            rotation_matrix=self.transform_matrix[:,1:4].reshape((3,3)),
            translate_vector=self.transform_matrix[:,0].reshape((1,3)))
          output_file.write(model_parser.transform_line())
      output_file.write('TER\nEND\n')
    return self.output_pdb

class UpdateCoords():
  def __init__(self, atom_line:str, file_format:str,
         rotation_matrix:np.ndarray, translate_vector:np.ndarray) -> None:
    self.atom_line = atom_line
    self.file_format = file_format.lower()

    assert rotation_matrix.shape == (3, 3)
    self.rotation_matrix = rotation_matrix

    assert translate_vector.shape == (1, 3)
    self.translation_vector = translate_vector

    match self.file_format:
      case '.pdb':
        # Legacy PDB format
        self.group_pdb = atom_line[0:6]    # Record type
        self.x = atom_line[30:38]      # X coordinate
        self.y = atom_line[38:46]      # Y coordinate
        self.z = atom_line[46:54]      # Z coordinate
      case '.cif':
        # mmCIF/PDBx files contain the same info as legacy PDB files
        # but only care about the order of each piece of info, rather
        # than the specific columns
        self.line_contents = [t for t in atom_line.split() if t != '']
        self.group_pdb = self.line_contents[0]       # _atom_site.group_PDB
        self.x = self.line_contents[10]          # _atom_site.Cartn_x
        self.y = self.line_contents[11]          # _atom_site.Cartn_y
        self.z = self.line_contents[12]          # _atom_site.Cartn_z
      case _:
        raise ValueError(f'Invalid file format: {self.file_format}. Expected .pdb or .cif')

  def _get_coords(self) -> np.ndarray:
    return np.array([float(self.x), float(self.y), float(self.z)])

  def _transform_coords(self) -> np.ndarray:
    return np.dot(self.rotation_matrix, self._get_coords()) + self.translation_vector

  def transform_line(self) -> str:
    transformed_coords = self._transform_coords()
    match self.file_format:
      case '.pdb':
        return self._pdb_format(transformed_coords)
      case '.cif':
        return self._cif_format(transformed_coords)
      case _:
        return ''

  def _pdb_format(self, transformed_coords:np.ndarray) -> str:
    return (
      f'{self.atom_line[:30]}'
      f'{transformed_coords[0][0]:8.3f}{transformed_coords[0][1]:8.3f}{transformed_coords[0][2]:8.3f}'
      f'{self.atom_line[54:]}'
    )

  def _cif_format(self, transformed_coords:np.ndarray) -> str:
    self.line_contents[10] = f'{transformed_coords[0][0]:8.3f}'
    self.line_contents[11] = f'{transformed_coords[0][1]:8.3f}'
    self.line_contents[12] = f'{transformed_coords[0][2]:8.3f}'
    return '\t'.join(self.line_contents)

class ColabAlign():
  def __init__(self, script_args) -> None:
    # Handle model list generation from user input
    self.model_list = []
    for user_input in script_args.input:
      if not user_input.is_file and not user_input.is_dir:
        continue
      if user_input.is_file and user_input.suffix in ('.pdb', '.cif'):
        self.model_list.append(user_input)
      elif user_input.is_dir():
        self.model_list.extend(
          file for file in user_input.iterdir() if file.is_file() and file.suffix in ('.pdb', '.cif')
        )
    self.model_list = sorted(self.model_list)

    # Handle user-defined clustering threshold
    self.thresholds = script_args.threshold

    # Handle input and output paths and create necessary directories
    print('Setting up file structure.')
    self.home_path = script_args.output
    self.models_path = self.home_path.joinpath('models')
    self.results_path = self.home_path.joinpath('results')

    for output_subdir in (
      self.home_path,
      self.models_path,
      self.results_path,
      self.results_path.joinpath('clusters'),
      ):
      Path.mkdir(output_subdir, exist_ok=True, parents=True)

    # Handle user-defined core count
    self.cores = script_args.cores
    if self.cores <= 0:
      warnings.warn(
        'User-defined core count cannot be lower than 1. Core count set to default (1).'
      )
      self.cores = 1
    if self.cores > cpu_count():
      warnings.warn(
        f'User-defined core count ({self.cores}) exceeds available cores ({cpu_count()}).'
        'Using maximum available cores instead.'
      )
      self.cores = cpu_count()

    self.usalign_path = script_args.path
    self.tmmatrix = None
    self.usalign_df = None

  def _reverse_transformation_matrix(self, transform_mx_forward):
    assert isinstance(transform_mx_forward, np.ndarray)
    assert transform_mx_forward.shape == (3, 4)

    translate_v_forward = transform_mx_forward[:, 0:1]
    rotation_mx_forward = transform_mx_forward[:, 1:4]

    rotation_mx_reverse = rotation_mx_forward.T
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

    transform_mx_forward = np.array([[d for d in line.split(' ') if d][1:] for line in decoded_results[4:7]]).astype(np.float64)
    transform_mx_reverse = self._reverse_transformation_matrix(transform_mx_forward.astype(np.float64))
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

  # Make phylogenetic tree from US-align results
  def pairwise_alignment(self):
    def _chunks(lst, n):
      '''Yield n chunks from the list.'''
      avg = len(lst) / n
      last = 0
      while last < len(lst):
        yield lst[int(last):int(last + avg)]
        last += avg

    print('Creating hardlinks for model files.')
    for model in self.model_list:
      hardlink_path = self.models_path.joinpath(model.name)
      if not hardlink_path.exists():
        hardlink_path.hardlink_to(model.resolve())

    # Calculate all possible combinations of models, then create discrete lists of comparisons
    # on a per-model basis. This enables us to run multiple, concurrent US-align instances and
    # avoiding doing any redundant calculations
    all_combos = combinations(self.model_list, 2)

    # Run multiple US-align jobs in parallel
    # This greatly speeds up the workflow
    print('Running US-align jobs.')
    process_results = []  # Collect dictionaries from each future

    with ProcessPoolExecutor(max_workers=self.cores) as executor:
      futures = [executor.submit(self._usalign_process, job) for job in _chunks(list(all_combos), self.cores)]
      for future in as_completed(futures):
        try:
          # Default dict set up to make a new empty subdict if the key does not exist
          process_hashmap = defaultdict(lambda: defaultdict(dict))

          for result in future.result():
            model1, model2, stdout, stderr = result

            # Parse the stdout from USalign
            tm_1, tm_2, transform_mx_forward, transform_mx_reverse = self._parse_usalign_stdout(stdout)

            # Update the hashmap with results
            process_hashmap.setdefault(model1.stem, {})[model2.stem] = {
              'tm_forward': tm_1,
              'transform_mx_forward': transform_mx_forward,
              'tm_reverse': tm_2,
              'transform_mx_reverse': transform_mx_reverse
            }

            # Handle any errors in stderr
            if stderr != b'' or stdout == b'':
              print(f"Error for models {model1, model2}: {stderr.decode(errors='ignore')}")

        except Exception as e:
          print(f"Error occurred: {e}")

        # Append this process's results
        process_results.append(process_hashmap)

    # Merge all results into a final dictionary
    usalign_results_map = {}
    for submap in process_results:
      for model1, subsubmap in submap.items():
        usalign_results_map.setdefault(model1, {}).update(subsubmap)

    # Export this to file?
    self.usalign_df = self._results_map_to_df(usalign_results_map)

    # Initialise a matrix with NaN values to be populated with US-align values
    # Probably not the most memory/time efficient method for this but I found it much easier
    # to understand how the data is handled.
    matrix_data = np.full((len(self.model_list), len(self.model_list)), np.nan)

    # Create a dictionary to map model names to indices
    model_index = {model.stem: idx for idx, model in enumerate(self.model_list)}

    # Populate the matrix with the highest US-align score for each pair
    # Chain identifiers are added after each .pdb or .cif but the first
    # chain in the file is the only one included in the alignment.
    # These identifiers are removed here to save confusion.
    for _, row in self.usalign_df.iterrows():
      if row['model1'].startswith('model1'):
        continue

      max_tm = max(row['tm_forward'], row['tm_reverse'])

      i, j = model_index[row['model1']], model_index[row['model2']]
      matrix_data[i, j] = max_tm
      matrix_data[j, i] = max_tm

    # Fill the diagonal with scores of 1
    # We don't need to align a model to itself because this will always return a TM-score of 1
    # Calculated US-align scores are given to 4 decimal places
    np.fill_diagonal(matrix_data, 1.0000)

    self.tmmatrix = pd.DataFrame(matrix_data, index=[model.stem for model in self.model_list], columns=[model.stem for model in self.model_list])
    self.tmmatrix.to_csv(self.results_path.joinpath('us-align_score_matrix.csv'), index=True)

    print('Generating phylogenetic tree.')
    # Invert US-align scores to make them suitable for distances on a phylogenetic tree
    # More similar pairs of models (i.e., higher TM-scores) have shorter distances to each other
    # Distance values are all rounded to 4 decimal places since all US-align scores are also to 4 d.p.
    # Clipping is also applied to ensure no values below 0 or above 1 are present in the distance matrix

    distances_df = 1.0000 - self.tmmatrix
    distances_df = distances_df.round(4)
    distances_df = distances_df.clip(lower=0.0000,upper=1.0000)

    # Convert the scores matrix to a lower triangle matrix
    # The lower and upper triangles of the matrix are identical and Bio.Phylo.TreeConstruction
    # requires a lower triangle matrix rather than the full matrix
    lower_tri_df = distances_df.where(np.tril(np.ones(distances_df.shape)).astype(bool))
    lower_tri_lists = [[value for value in row if not isnan(value)] for row in lower_tri_df.values.tolist()]

    # Generate phylogenetic tree using the UPGMA clustering method
    # We can safely ignore the Molecular Clock hypothesis because we are not
    # deriving evolutionary relationships between proteins
    tm_matrix = DistanceMatrix(names=[model.stem for model in self.model_list], matrix=lower_tri_lists)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(tm_matrix)

    # Draw the tree in ASCII for quick validation
    # print('\n')
    # Phylo.draw_ascii(tree)
    Phylo.write(tree, self.results_path.joinpath('colabalign.tree'), 'newick')

  def tree_clustering(self):
    print('Calculating clusters from structural tree.')

    completed_processes = []

    for thresh in self.thresholds:
      cmd = [
        'TreeCluster.py',
           '-i', self.results_path.joinpath('colabalign.tree').as_posix(),
           '-o', self.results_path.joinpath(f'clusters/{thresh:.2f}.tsv').as_posix(),
          '-t', f'{thresh:.2f}'
          ]
      with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
        _, _ = process.communicate()
        completed_processes.append(process)

    for process in completed_processes:
      process.wait()

    # Make lists of models for each cluster
    clusters = {}

    for thresh in self.thresholds:
      with open(self.results_path.joinpath(f'clusters/{thresh:.2f}.tsv'), 'r', encoding='UTF8') as cluster_file:
        for line in cluster_file:
          if line.startswith('SequenceName'):
            continue
          line_split = line.split('\t')
          cluster_num = int(line_split[1].strip('\n'))
          if cluster_num not in clusters:
            clusters[cluster_num] = []
          clusters[cluster_num].append(line_split[0].strip())

      print(f'Clustering Threshold: {thresh:.2f}')
      for cluster_num, model_list in clusters.items():
        print(f'Aligning cluster: {cluster_num}')
        cluster_output_dir = self.results_path.joinpath(f'clusters/{thresh:.2f}/{cluster_num}')
        Path.mkdir(cluster_output_dir, exist_ok=True, parents=True)

        if cluster_num == -1:
          for model in model_list:
            copy(
              src=self.models_path.joinpath(f'{model}.pdb'),
              dst=cluster_output_dir.joinpath(f'{model}.pdb')
              )
          continue

        sub_matrix = self.tmmatrix.loc[list(model_list), list(model_list)]

        reference_model_name = sub_matrix.max().idxmax()

        with (open(self.models_path.joinpath(f'{reference_model_name}.pdb'), 'r', encoding='UTF8') as input_ref,
              open(cluster_output_dir.joinpath(f'{reference_model_name}.pdb'), 'w', encoding='UTF8') as output_ref):
            output_ref.write(f'REMARK 900\nREMARK 900 RELATED ENTRIES\nREMARK 900 REFERENCE MODEL FOR CLUSTER {cluster_num}\n')
            output_ref.write(input_ref.read())

        cluster_df = self.usalign_df[(self.usalign_df['model1'] == reference_model_name) | (self.usalign_df['model2'] == reference_model_name)]

        for model in model_list:
          if model == reference_model_name:
            continue

          mx = cluster_df.apply(
            self._find_matrix_in_array,
            axis=1,
            args=(model,)
            ).dropna().iloc[0]

          if mx is not None:
            aligner = StructureAligner(input_pdb=self.models_path.joinpath(f'{model}.pdb'),
                        output_pdb=cluster_output_dir.joinpath(f'{model}.pdb'),
                        transform_matrix=mx)
            aligner.transform_coords()

# For running locally
def local():
  local_instance = ColabAlign(script_args())
  local_instance.pairwise_alignment()
  local_instance.tree_clustering()

# For running on Google Colab
def colab(args):
  # The colab notebook sets the user arguments with sys.argv
  colab_instance = ColabAlign(args)
  colab_instance.pairwise_alignment()
  colab_instance.tree_clustering()

if __name__ == '__main__':
  local()
