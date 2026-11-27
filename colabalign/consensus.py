'''Progressive profile-merge MSA construction and consensus/representative selection.'''

from math import inf
from pathlib import Path
from shutil import copyfile
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from numba import njit as _numba_njit
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from tqdm.auto import tqdm as tqdm_auto

from .common import extract_sequence, process_context

NEG_INF = -1e9


def njit(function):
    '''Compile with Numba, caching the machine code on disk when possible.

    Numba refuses to cache when neither the package folder nor the user cache
    folder is writable (e.g. a read-only container); compile per process then.'''
    try:
        return _numba_njit(cache=True, nogil=True)(function)
    except RuntimeError:
        return _numba_njit(nogil=True)(function)

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
}




# Progressive profile-merge MSA (per cluster)
def _column_luts(profile):
    '''For every structure, an array mapping residue index -> profile column (-1 if absent).

    Structures are listed in order of first appearance in the profile, which
    sets the (floating-point) order in which pair scores are accumulated.'''
    positions = {}
    for col_idx, col in enumerate(profile):
        for sid, res in col.items():
            entry = positions.get(sid)
            if entry is None:
                positions[sid] = entry = ([], [])
            entry[0].append(res)
            entry[1].append(col_idx)
    luts = {}
    for sid, (residues, columns) in positions.items():
        lut = np.full(max(residues) + 1, -1, dtype=np.int32)
        lut[residues] = columns  # each residue occurs once in a profile
        luts[sid] = lut
    return luts


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
    luts_a, luts_b = _column_luts(profile_a), _column_luts(profile_b)

    for sid_a, lut_a in luts_a.items():
        for sid_b, lut_b in luts_b.items():
            res_a, res_b, weight = _resolve_pair_entry(sid_a, sid_b, library)
            if res_a is None or len(res_a) == 0:
                continue
            _accumulate_pair_score(score, lut_a, lut_b, res_a, res_b, weight)

    return score


@njit
def _fill_tables(score, gap_open, gap_extend, neg_inf):
    '''Fill the three Gotoh tables (compiled with Numba).

    Each max keeps the first argument on a tie, as Python's max() does.'''
    n, m = score.shape
    match = np.full((n + 1, m + 1), neg_inf)
    gap_x = np.full((n + 1, m + 1), neg_inf)
    gap_y = np.full((n + 1, m + 1), neg_inf)
    match[0, 0] = 0.0
    for i in range(1, n + 1):
        gap_x[i, 0] = -gap_open - (i - 1) * gap_extend
    for j in range(1, m + 1):
        gap_y[0, j] = -gap_open - (j - 1) * gap_extend
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            best = match[i-1, j-1]
            if gap_x[i-1, j-1] > best:
                best = gap_x[i-1, j-1]
            if gap_y[i-1, j-1] > best:
                best = gap_y[i-1, j-1]
            match[i, j] = best + np.float64(score[i-1, j-1])
            opened = match[i-1, j] - gap_open
            extended = gap_x[i-1, j] - gap_extend
            gap_x[i, j] = extended if extended > opened else opened
            opened = match[i, j-1] - gap_open
            extended = gap_y[i, j-1] - gap_extend
            gap_y[i, j] = extended if extended > opened else opened
    return match, gap_x, gap_y


def gotoh_align(score, gap_open=6.0, gap_extend=1.0):
    '''Run Gotoh affine-gap alignment on a score matrix and return the traceback path.

    M[i,j] = max(M, X, Y)[i-1,j-1] + S[i-1,j-1]
    X[i,j] = max(M[i-1,j] - gap_open, X[i-1,j] - gap_extend)
    Y[i,j] = max(M[i,j-1] - gap_open, Y[i,j-1] - gap_extend)'''
    match, gap_x, gap_y = _fill_tables(np.ascontiguousarray(score, dtype=np.float32),
                                       float(gap_open), float(gap_extend), float(NEG_INF))
    rows, cols = _traceback(match, gap_x, gap_y)
    return [(None if i < 0 else int(i), None if j < 0 else int(j)) for i, j in zip(rows, cols)]


@njit
def _best_state(match, gap_x, gap_y, i, j):
    '''Index of the largest of M, X, Y at (i, j); the first wins a tie, as np.argmax.'''
    state, best = 0, match[i, j]
    if gap_x[i, j] > best:
        state, best = 1, gap_x[i, j]
    if gap_y[i, j] > best:
        state = 2
    return state


@njit
def _traceback(match, gap_x, gap_y):
    '''Trace the optimal path back from the bottom-right corner (compiled with Numba).

    Returns the aligned row and column indices in path order, with -1 for a gap.'''
    i, j = match.shape[0] - 1, match.shape[1] - 1
    rows = np.empty(i + j, dtype=np.int64)
    cols = np.empty(i + j, dtype=np.int64)
    k = 0
    state = _best_state(match, gap_x, gap_y, i, j)
    while i > 0 or j > 0:
        if state == 0 and i > 0 and j > 0:
            rows[k], cols[k] = i - 1, j - 1
            i, j = i - 1, j - 1
            state = 0 if (i == 0 and j == 0) else _best_state(match, gap_x, gap_y, i, j)
        elif state == 1 and i > 0:
            rows[k], cols[k] = i - 1, -1
            i -= 1
            state = 0 if match[i, j] >= gap_x[i, j] else 1
        elif j > 0:
            rows[k], cols[k] = -1, j - 1
            j -= 1
            state = 0 if match[i, j] >= gap_y[i, j] else 2
        else:
            break
        k += 1
    return rows[:k][::-1], cols[:k][::-1]


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

def extract_leaf_profile(sid, profile):
    '''Extract a single structure's columns out of a merged profile as its own profile.'''
    return [{sid: col[sid]} for col in profile if sid in col]


def remove_structure(sid, profile):
    '''Remove a structure from a merged profile, dropping any now-empty columns.'''
    remaining = [{k: v for k, v in col.items() if k != sid} for col in profile]
    return [col for col in remaining if col]


def structure_score(sid, luts, library, others):
    '''Weights of the library pairs between sid and each other structure that the
    profile places in the same column. Summed over all sids, every pair is counted
    twice, so the sum over all structures is twice the profile's total score.'''
    lut_a = luts[sid]
    total = 0.0
    for other in others:
        res_a, res_b, weight = _resolve_pair_entry(sid, other, library)
        if res_a is None or len(res_a) == 0:
            continue
        lut_b = luts[other]
        valid = (res_a < len(lut_a)) & (res_b < len(lut_b))
        col_a, col_b = lut_a[res_a[valid]], lut_b[res_b[valid]]
        same = (col_a >= 0) & (col_a == col_b)
        total += float(weight[valid][same].sum(dtype=np.float64))
    return total


def refine_msa(profile, library, max_passes=3):
    '''Improve an MSA via leave-one-out re-alignment passes until no gain is found.

    Re-aligning one structure keeps every column of the remaining profile intact,
    so only the pairs involving that structure can change their score. Comparing
    those pairs alone gives the same accept/reject decision as re-scoring the
    whole candidate and current profiles, in O(|C|) rather than O(|C|^2) pair
    look-ups per candidate.'''
    struct_ids = sorted({sid for col in profile for sid in col})
    if len(struct_ids) < 3:
        return profile  # nothing to refine with fewer than 3 structures

    current = profile
    current_luts = _column_luts(current)

    for _ in range(max_passes):
        improved = False
        for sid in struct_ids:
            others = [other for other in struct_ids if other != sid]
            leaf = extract_leaf_profile(sid, current)
            rest = remove_structure(sid, current)
            candidate = merge_profiles(rest, leaf, library)
            candidate_luts = _column_luts(candidate)
            if (structure_score(sid, candidate_luts, library, others)
                    > structure_score(sid, current_luts, library, others)):
                current, current_luts = candidate, candidate_luts
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




def build_cluster_alignment(sequences, library, tmmatrix, refine=True):
    '''Build (and optionally refine) the aligned sequences for one cluster.'''
    members = list(sequences)
    if len(members) == 1:
        return sequences
    root_profile = build_cluster_msa(members, sequences, library, tmmatrix)
    if refine:
        root_profile = refine_msa(root_profile, library)
    return profile_to_aligned_sequences(root_profile, sequences)


def _write_cluster_outputs(pdb_dir, cluster_num, threshold, aligned, consensus):
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


def _copy_cluster_outputs(source_dir, threshold, cluster_num, target_dir, target_threshold,
                          target_num):  # pylint: disable=too-many-arguments,too-many-positional-arguments
    '''Copy an identical cluster's alignment and consensus under another threshold's names.'''
    for suffix in ('.afasta', '_consensus.afasta'):
        source = _cluster_file(source_dir, threshold, cluster_num, suffix)
        if source.exists():
            copyfile(source, _cluster_file(target_dir, target_threshold, target_num, suffix))


def _cluster_file(pdb_dir, threshold, cluster_num, suffix):
    '''Path of a per-cluster output file, e.g. threshold_0.20_cluster_3.afasta.'''
    return Path(pdb_dir).joinpath(f'threshold_{threshold}_cluster_{cluster_num}{suffix}')


# pylint: disable=too-many-arguments,too-many-positional-arguments
def process_cluster(pdb_dir: Path, cluster_num, threshold, pair_store, sub_tmmatrix,
                    refine=True, member_files=None):
    '''
    Build the MSA/consensus for one cluster and return its representative.
    Runs in a worker process. It receives only this cluster's slice of the
    TM-score matrix, and loads only this cluster's pairs from the PairStore.
    Sequences come from the superposed PDBs in pdb_dir, or from member_files
    (the prepared models) when superposition was skipped.
    '''
    pdb_files = sorted(member_files, key=lambda p: p.stem) if member_files \
        else sorted(pdb_dir.glob('*.pdb'))
    if not pdb_files:
        return None, cluster_num, threshold

    sequences = {p.stem: extract_sequence(p) for p in pdb_files}
    library = pair_store.sub_library(list(sequences))
    aligned = build_cluster_alignment(sequences, library, sub_tmmatrix, refine)
    consensus = build_consensus(aligned, threshold=0.6)
    _write_cluster_outputs(pdb_dir, cluster_num, threshold, aligned, consensus)

    representative = find_nearest_sequence(aligned, consensus)
    return representative, cluster_num, threshold


class ClusterConsensus:
    '''Builds per-cluster MSAs, consensus sequences, and selects representatives.'''

    # pylint: disable=too-many-arguments,too-many-positional-arguments
    def __init__(self, pair_store, tmmatrix, threads: int = 1, refine: bool = True,
                 models_by_name=None) -> None:
        '''Store the PairStore, TM-score matrix, and run parameters.

        models_by_name maps structure names to prepared model files. When given,
        sequences are read from these rather than from the superposed copies.'''
        self.pair_store = pair_store
        self.tmmatrix = tmmatrix
        self.threads = threads
        self.refine = refine
        self.models_by_name = models_by_name

    def _collect_jobs(self, results_path: Path, clusters_by_threshold: dict):
        '''Build the list of (threshold, cluster_num, pdb_dir, member names) jobs.'''
        jobs = []
        for threshold_str, clusters in clusters_by_threshold.items():
            for cluster_num, member_names in clusters.items():
                if cluster_num == -1:
                    continue
                pdb_dir = results_path.joinpath(f'clusters/{threshold_str}/{cluster_num}')
                jobs.append((threshold_str, cluster_num, pdb_dir, member_names))
        return jobs

    def _run_jobs(self, jobs):
        '''Run cluster-processing jobs in a shared process pool and collect results.

        Nearby thresholds share most of their clusters. Each distinct member set is
        aligned once; its alignment files are copied to the other thresholds.'''
        representatives_by_threshold = {threshold_str: {} for threshold_str, *_ in jobs}
        groups = {}
        for job in jobs:
            groups.setdefault(frozenset(job[3]), []).append(job)

        with ProcessPoolExecutor(max_workers=self.threads, mp_context=process_context()) as executor:
            total_bar = tqdm_auto(
                total=len(groups), desc='Total clusters',
                unit=' reps chosen', dynamic_ncols=True, leave=True,
            )
            futures = {}
            # Largest clusters first: one big cluster can take longer than all the
            # others together, so starting it last would leave workers idle.
            for member_set, group in sorted(groups.items(),
                                            key=lambda item: (-len(item[0]), sorted(item[0]))):
                threshold_str, cluster_num, pdb_dir, _ = group[0]
                members = [m for m in sorted(member_set) if m in self.tmmatrix.index]
                sub_tmmatrix = self.tmmatrix.loc[members, members]
                member_files = ([self.models_by_name[m] for m in sorted(member_set)]
                                if self.models_by_name else None)
                pdb_dir.mkdir(parents=True, exist_ok=True)
                futures[executor.submit(
                    process_cluster, pdb_dir, cluster_num, threshold_str,
                    self.pair_store, sub_tmmatrix, self.refine, member_files,
                )] = group
            for future in as_completed(futures):
                representative, cluster_num, threshold_str = future.result()
                total_bar.update(1)
                if representative is None:
                    continue
                group = futures[future]
                source_dir = group[0][2]
                for other_threshold, other_num, other_dir, _ in group:
                    if (other_threshold, other_num) != (threshold_str, cluster_num):
                        _copy_cluster_outputs(source_dir, threshold_str, cluster_num,
                                              other_dir, other_threshold, other_num)
                    representatives_by_threshold[other_threshold][int(other_num)] = representative
        total_bar.close()

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
        '''Run a single shared process pool across every threshold's clusters, write the
        results, and return dict[threshold_str -> dict[cluster_num -> representative]].'''
        jobs = self._collect_jobs(results_path, clusters_by_threshold)
        if not jobs:
            return {threshold_str: {} for threshold_str in clusters_by_threshold}
        representatives_by_threshold = self._run_jobs(jobs)
        self._write_representatives(results_path, representatives_by_threshold)
        for threshold_str in clusters_by_threshold:
            representatives_by_threshold.setdefault(threshold_str, {})
        return representatives_by_threshold
