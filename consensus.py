'''Progressive profile-merge MSA construction and consensus/representative selection.'''

from math import inf
from pathlib import Path
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from tqdm.auto import tqdm as tqdm_auto

from common import extract_sequence

NEG_INF = -1e9

MVIEW_LUT = {
    'o': {'S', 'T'},
    'l': {'I', 'L', 'V'},
    'a': {'F', 'H', 'W', 'Y'},
    'c': {'D', 'E', 'H', 'K', 'R'},
    'h': {'A', 'C', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'R', 'T', 'V', 'W', 'Y'},
    '-': {'D', 'E'},
    'p': {'C', 'D', 'E', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'},
    '+': {'H', 'K', 'R'},
    's': {'A', 'C', 'D', 'G', 'N', 'P', 'S', 'T', 'V'},
    'u': {'A', 'G', 'S'},
    't': {'A', 'C', 'D', 'E', 'G', 'H', 'K', 'N', 'Q', 'R', 'S', 'T'},
    '*': {'*'},
}




# Progressive profile-merge MSA (per cluster)
def _residue_to_col_map(profile):
    '''Map each structure ID to a {residue_index: column_index} dict for a profile.'''
    maps = {}
    for col_idx, col in enumerate(profile):
        for sid, res in col.items():
            maps.setdefault(sid, {})[res] = col_idx
    return maps


def _build_residue_lut(residue_to_col):
    '''Build a dense lookup array mapping residue index -> profile column index.'''
    max_res = max(residue_to_col.keys(), default=-1) + 1
    lut = np.full(max_res, -1, dtype=np.int32)
    for residue, col in residue_to_col.items():
        lut[residue] = col
    return lut


def _resolve_pair_entry(sid_a, sid_b, library):
    '''Look up the alignment entry for a structure pair, handling key order.'''
    key = (sid_a, sid_b) if (sid_a, sid_b) in library else (sid_b, sid_a)
    entry = library.get(key)
    if entry is None:
        return None, None, None
    flip = key != (sid_a, sid_b)
    res_a = entry['res_b'] if flip else entry['res_a']
    res_b = entry['res_a'] if flip else entry['res_b']
    return res_a, res_b, entry['weight']


def _accumulate_pair_score(score, lut_a, lut_b, res_a, res_b, weight):
    '''Add pairwise alignment weights into the score matrix for one structure pair.'''
    valid = (res_a < len(lut_a)) & (res_b < len(lut_b))
    col_i = lut_a[res_a[valid]]
    col_j = lut_b[res_b[valid]]
    valid2 = (col_i >= 0) & (col_j >= 0)
    if np.any(valid2):
        np.add.at(score, (col_i[valid2], col_j[valid2]), weight[valid][valid2])


def build_score_matrix_sparse(profile_a, profile_b, library):
    '''Build a dense column-vs-column score matrix between two profiles.'''
    score = np.zeros((len(profile_a), len(profile_b)), dtype=np.float32)
    map_a, map_b = _residue_to_col_map(profile_a), _residue_to_col_map(profile_b)

    for sid_a, residues_a in map_a.items():
        lut_a = _build_residue_lut(residues_a)
        for sid_b, residues_b in map_b.items():
            res_a, res_b, weight = _resolve_pair_entry(sid_a, sid_b, library)
            if res_a is None or len(res_a) == 0:
                continue
            lut_b = _build_residue_lut(residues_b)
            _accumulate_pair_score(score, lut_a, lut_b, res_a, res_b, weight)

    return score


def gotoh_align(score, gap_open=6.0, gap_extend=1.0):
    '''Run Gotoh affine-gap alignment on a score matrix and return the traceback path.'''
    n, m = score.shape
    match = np.full((n + 1, m + 1), NEG_INF)
    gap_x = np.full((n + 1, m + 1), NEG_INF)
    gap_y = np.full((n + 1, m + 1), NEG_INF)
    match[0, 0] = 0.0
    for i in range(1, n + 1):
        gap_x[i, 0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        gap_y[0, j] = -gap_open - (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match[i, j] = max(match[i-1, j-1], gap_x[i-1, j-1], gap_y[i-1, j-1]) + score[i-1, j-1]
            gap_x[i, j] = max(match[i-1, j] - gap_open, gap_x[i-1, j] - gap_extend)
            gap_y[i, j] = max(match[i, j-1] - gap_open, gap_y[i, j-1] - gap_extend)

    i, j = n, m
    state = int(np.argmax([match[i, j], gap_x[i, j], gap_y[i, j]]))
    path = []
    while i > 0 or j > 0:
        if state == 0 and i > 0 and j > 0:
            path.append((i - 1, j - 1))
            i, j = i - 1, j - 1
            state = 0 if (i, j) == (0, 0) else int(
                np.argmax([match[i, j], gap_x[i, j], gap_y[i, j]])
                )
        elif state == 1 and i > 0:
            path.append((i - 1, None))
            i -= 1
            state = 0 if match[i, j] >= gap_x[i, j] else 1
        elif j > 0:
            path.append((None, j - 1))
            j -= 1
            state = 0 if match[i, j] >= gap_y[i, j] else 2
        else:
            break
    path.reverse()
    return path


def merge_profiles(profile_a, profile_b, library, gap_open=6.0, gap_extend=1.0):
    '''Align and merge two profiles into a single combined profile.'''
    score = build_score_matrix_sparse(profile_a, profile_b, library)
    path = gotoh_align(score, gap_open, gap_extend)
    merged = []
    for i, j in path:
        col = {}
        if i is not None:
            col.update(profile_a[i])
        if j is not None:
            col.update(profile_b[j])
        merged.append(col)
    return merged


def build_cluster_guide_tree(members, tmmatrix):
    '''Build a hierarchical-clustering guide tree from pairwise TM-scores.'''
    n = len(members)
    dist = np.zeros((n, n), dtype=np.float64)
    for i, a in enumerate(members):
        for j, b in enumerate(members):
            if i < j:
                d = 1.0 - float(tmmatrix.loc[a, b])
                dist[i, j] = dist[j, i] = max(d, 0.0)
    return linkage(squareform(dist, checks=False), method='average')


def build_cluster_msa(members, sequences, library, tmmatrix):
    '''Progressively build a multiple sequence alignment for a cluster of structures.'''
    n = len(members)
    profiles = {i: [{members[i]: r} for r in range(len(sequences[members[i]]))] for i in range(n)}

    if n == 1:
        return profiles[0]
    if n == 2:
        return merge_profiles(profiles[0], profiles[1], library)

    guide_tree = build_cluster_guide_tree(members, tmmatrix)
    next_id = n
    for row in guide_tree:
        left, right = int(row[0]), int(row[1])
        profiles[next_id] = merge_profiles(profiles[left], profiles[right], library)
        del profiles[left], profiles[right]
        next_id += 1
    return profiles[next_id - 1]


def profile_to_aligned_sequences(root_profile, sequences):
    '''Convert a merged profile back into per-structure aligned sequence strings.'''
    aligned = {sid: [] for sid in sequences}
    for col in root_profile:
        for sid in sequences:
            aligned[sid].append(sequences[sid][col[sid]] if sid in col else '-')
    return {sid: ''.join(chars) for sid, chars in aligned.items()}

def profile_score_sum(profile, library):
    '''Sum pairwise alignment weights realized by the current column arrangement.'''
    struct_ids = sorted({sid for col in profile for sid in col})
    total = 0.0
    for a_idx, sid_a in enumerate(struct_ids):
        for sid_b in struct_ids[a_idx + 1:]:
            res_a, res_b, weight = _resolve_pair_entry(sid_a, sid_b, library)
            if res_a is None:
                continue
            aligned_pairs = {
                (ra, rb): wi for ra, rb, wi in zip(res_a.tolist(), res_b.tolist(), weight.tolist())
            }
            for col in profile:
                if sid_a in col and sid_b in col:
                    total += aligned_pairs.get((col[sid_a], col[sid_b]), 0.0)
    return total


def extract_leaf_profile(sid, profile):
    '''Extract a single structure's columns out of a merged profile as its own profile.'''
    return [{sid: col[sid]} for col in profile if sid in col]


def remove_structure(sid, profile):
    '''Remove a structure from a merged profile, dropping any now-empty columns.'''
    remaining = [{k: v for k, v in col.items() if k != sid} for col in profile]
    return [col for col in remaining if col]


def refine_msa(profile, library, max_passes=3):
    '''Improve an MSA via leave-one-out re-alignment passes until no gain is found.'''
    struct_ids = sorted({sid for col in profile for sid in col})
    if len(struct_ids) < 3:
        return profile  # nothing to refine with fewer than 3 structures

    current = profile
    current_score = profile_score_sum(current, library)

    for _ in range(max_passes):
        improved = False
        for sid in struct_ids:
            leaf = extract_leaf_profile(sid, current)
            rest = remove_structure(sid, current)
            candidate = merge_profiles(rest, leaf, library)
            candidate_score = profile_score_sum(candidate, library)
            if candidate_score > current_score:
                current, current_score = candidate, candidate_score
                improved = True
        if not improved:
            break
    return current




# MView-style consensus + representative selection

def consensus_column(chars, threshold=0.6):
    '''Pick a consensus symbol (residue or MView-style group code) for one column.'''
    non_gap = [c for c in chars if c != '-']
    if not non_gap:
        return '.'
    n = len(non_gap)
    top_char, top_count = Counter(non_gap).most_common(1)[0]
    if top_count / n >= threshold:
        return top_char
    for symbol, group in MVIEW_LUT.items():
        if symbol == '*':
            continue
        if sum(1 for c in non_gap if c in group) / n >= threshold:
            return symbol
    return '.'


def build_consensus(aligned_sequences, threshold=0.6):
    '''Build a full consensus sequence string from a set of aligned sequences.'''
    seqs = list(aligned_sequences.values())
    length = len(seqs[0])
    return ''.join(consensus_column([s[i] for s in seqs], threshold) for i in range(length))


def _is_match(c1, c2):
    '''Check whether an aligned residue matches a consensus symbol or group code.'''
    if c1 == '.' and c2 == '-':
        return True
    if c1 == '.':
        return False
    if c1 == c2:
        return True
    return c2 in MVIEW_LUT.get(c1, set())


def find_nearest_sequence(aligned_sequences, consensus):
    '''Find the aligned sequence with the highest match score against the consensus.'''
    best = {'record_name': '', 'sequence': '', 'score': -inf, 'max_score': len(consensus)}
    for name, seq in aligned_sequences.items():
        score = sum(
            1 if _is_match(cons_char, seq_char) else -1
            for cons_char, seq_char in zip(consensus, seq)
        )
        if score > best['score']:
            best = {
                'record_name': name, 'sequence': seq,
                'score': score, 'max_score': len(consensus),
            }
    return best


def write_afasta(aligned_sequences, path, width=80):
    '''Write aligned sequences to a FASTA file, wrapped at the given line width.'''
    with open(path, 'w', encoding='utf-8') as f:
        for name, seq in aligned_sequences.items():
            f.write(f'>{name}\n')
            for i in range(0, len(seq), width):
                f.write(f'{seq[i:i+width]}\n')




class ClusterConsensus:
    '''Builds per-cluster MSAs, consensus sequences, and selects representatives.'''

    def __init__(self, library: dict, tmmatrix, threads: int = 1, refine: bool = True) -> None:
        '''Store shared alignment library, TM-score matrix, and run parameters.'''
        self.library = library
        self.tmmatrix = tmmatrix
        self.threads = threads
        self.refine = refine

    def _build_cluster_alignment(self, sequences):
        '''Build (and optionally refine) the aligned sequences for a cluster.'''
        members = list(sequences)
        if len(members) == 1:
            return sequences
        root_profile = build_cluster_msa(members, sequences, self.library, self.tmmatrix)
        if self.refine:
            root_profile = refine_msa(root_profile, self.library)
        return profile_to_aligned_sequences(root_profile, sequences)

    def _write_cluster_outputs(self, pdb_dir, cluster_num, threshold, aligned, consensus):
        '''Write the alignment and consensus FASTA files for one cluster.'''
        afasta_path = pdb_dir.joinpath(f'threshold_{threshold}_cluster_{cluster_num}.afasta')
        write_afasta(aligned, afasta_path)

        consensus_path = pdb_dir.joinpath(
            f'threshold_{threshold}_cluster_{cluster_num}_consensus.afasta'
        )
        with open(consensus_path, 'w', encoding='utf-8') as f:
            f.write('>consensus\n')
            for i in range(0, len(consensus), 80):
                f.write(f'{consensus[i:i+80]}\n')

    def _process_cluster(self, pdb_dir: Path, cluster_num, threshold):
        '''Build the MSA/consensus for one cluster and return its representative.'''
        pdb_files = sorted(pdb_dir.glob('*.pdb'))
        if not pdb_files:
            return None, cluster_num, threshold

        sequences = {p.stem: extract_sequence(p) for p in pdb_files}
        aligned = self._build_cluster_alignment(sequences)
        consensus = build_consensus(aligned, threshold=0.6)
        self._write_cluster_outputs(pdb_dir, cluster_num, threshold, aligned, consensus)

        representative = find_nearest_sequence(aligned, consensus)
        return representative, cluster_num, threshold

    def _collect_jobs(self, results_path: Path, clusters_by_threshold: dict):
        '''Build the list of (threshold, cluster_num, pdb_dir) jobs to process.'''
        jobs = []
        for threshold_str, clusters in clusters_by_threshold.items():
            for cluster_num in clusters:
                if cluster_num == -1:
                    continue
                pdb_dir = results_path.joinpath(f'clusters/{threshold_str}/{cluster_num}')
                jobs.append((threshold_str, cluster_num, pdb_dir))
        return jobs

    def _run_jobs(self, jobs):
        '''Run cluster-processing jobs in a shared process pool and collect results.'''
        representatives_by_threshold = {threshold_str: {} for threshold_str, _, _ in jobs}

        with ProcessPoolExecutor(max_workers=self.threads) as executor:
            total_bar = tqdm_auto(
                total=len(jobs), desc='Total clusters',
                unit=' reps chosen', dynamic_ncols=True, leave=True,
            )
            futures = [
                executor.submit(self._process_cluster, pdb_dir, cluster_num, threshold_str)
                for threshold_str, cluster_num, pdb_dir in jobs
            ]
            for fut in futures:
                fut.add_done_callback(lambda _: total_bar.update(1))

            for future in as_completed(futures):
                representative, cluster_num, threshold_str = future.result()
                if representative is not None:
                    representatives_by_threshold[threshold_str][int(cluster_num)] = representative

        return representatives_by_threshold

    def _write_representatives(self, results_path: Path, representatives_by_threshold: dict):
        '''Write the final per-threshold cluster-representative FASTA files.'''
        for threshold_str, representatives in representatives_by_threshold.items():
            if not representatives:
                continue
            out_path = results_path.joinpath(
                f'clusters/{threshold_str}/threshold_{threshold_str}_cluster_representatives.fasta'
            )
            ordered = dict(sorted(representatives.items(), key=lambda x: int(x[0])))
            with open(out_path, 'w', encoding='UTF8') as f:
                for num, entry in ordered.items():
                    f.write(
                        f'> Cluster {num}: {entry["record_name"]} '
                        f'(alignment score {entry["score"]}/{entry["max_score"]})\n'
                        f'{entry["sequence"].replace("-","")}\n'
                    )

    def get_representatives_all_thresholds(self, results_path: Path, clusters_by_threshold: dict):
        '''Run a single shared process pool across every threshold's clusters and write results.'''
        jobs = self._collect_jobs(results_path, clusters_by_threshold)
        if not jobs:
            return
        representatives_by_threshold = self._run_jobs(jobs)
        self._write_representatives(results_path, representatives_by_threshold)
