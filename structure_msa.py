#!/usr/bin/env python3
'''
Structure-informed MSA via USalign pairwise alignments + progressive profile merging.

Usage:
    python structure_msa.py --pdb-dir /path/to/pdbs [--cache-dir ./usalign_cache] [--workers 4]

Expects one PDB file per structure in --pdb-dir (single chain each).
'''

import argparse
import hashlib
import pickle
import re
import subprocess
from pathlib import Path
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

from common import extract_sequence

NEG_INF = -1e9


# USalign run + parse + disk cache
def _file_hash(path):
    '''Return a short sha256 hex digest of a file's contents, for cache keys.'''
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def _cache_key(path_a, path_b):
    '''Build a cache filename from the hashes of two input structure files.'''
    return f'{_file_hash(path_a)}_{_file_hash(path_b)}.pkl'


def _find_alignment_block(lines):
    '''Locate the (seq_a, match_line, seq_b) block within USalign stdout.'''
    valid_chars = set('ABCDEFGHIKLMNPQRSTVWYXZUO-')
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and set(stripped) <= valid_chars and len(stripped) > 3:
            return i
    return None


def _parse_do_distances(lines):
    '''Extract per-residue-pair CA distances from USalign's -do output lines.'''
    do_pattern = re.compile(
        r'^\s*CA\s+\w{3}\s+\S+\s+\d+\s+CA\s+\w{3}\s+\S+\s+\d+\s+([\d.]+)\s*$'
    )
    return [float(m.group(1)) for line in lines if (m := do_pattern.match(line))]


def _build_residue_pairs(seq_a, match_line, seq_b):
    '''Walk the alignment block to build aligned residue-index pairs and weights.'''
    res_a_idx = res_b_idx = 0
    pairs_a, pairs_b, weights = [], [], []
    for ca, m, cb in zip(seq_a, match_line, seq_b):
        a_gap, b_gap = ca == '-', cb == '-'
        if not a_gap and not b_gap and m in (':', '.'):
            pairs_a.append(res_a_idx)
            pairs_b.append(res_b_idx)
            weights.append(1.0 if m == ':' else 0.5)  # fallback binary/graded weight
        if not a_gap:
            res_a_idx += 1
        if not b_gap:
            res_b_idx += 1
    return pairs_a, pairs_b, np.array(weights, dtype=np.float32)


def _weight_by_distance(weights, do_distances):
    '''Recompute weights from CA distances (TM-score-style falloff) when available.'''
    if len(do_distances) != len(weights):
        # silently keep the binary/graded fallback rather than failing the whole pair
        return weights
    d = np.array(do_distances, dtype=np.float32)
    d0 = 0.5 if len(d) < 15 else max(1.24 * (len(d) - 15) ** (1 / 3) - 1.8, 0.5)
    return 1.0 / (1.0 + (d / d0) ** 2)


def parse_usalign_output(stdout_text):
    '''Parse TM-scores, RMSD, and the aligned residue pairing out of USalign stdout.'''
    tm_matches = re.findall(r'TM-score=\s*([\d.]+)', stdout_text)
    if len(tm_matches) < 2:
        raise ValueError('could not parse TM-scores from USalign output:\n' + stdout_text[:1000])
    tm1, tm2 = float(tm_matches[0]), float(tm_matches[1])

    rmsd_match = re.search(r'RMSD=\s*([\d.]+)', stdout_text)
    rmsd = float(rmsd_match.group(1)) if rmsd_match else None

    lines = stdout_text.splitlines()
    aln_start = _find_alignment_block(lines)
    if aln_start is None or aln_start + 2 >= len(lines):
        raise ValueError('could not locate alignment block in USalign output:\n' + stdout_text[:1000])

    seq_a, match_line, seq_b = lines[aln_start], lines[aln_start + 1], lines[aln_start + 2]
    do_distances = _parse_do_distances(lines)

    pairs_a, pairs_b, weights = _build_residue_pairs(seq_a, match_line, seq_b)
    weights = _weight_by_distance(weights, do_distances)

    return {
        'tm1': tm1, 'tm2': tm2, 'rmsd': rmsd,
        'res_a': np.array(pairs_a, dtype=np.int32),
        'res_b': np.array(pairs_b, dtype=np.int32),
        'weight': weights,
    }


def run_usalign_pair_cached(id_a, path_a, id_b, path_b, workdir):
    '''Run (or load from cache) a single pairwise USalign alignment.'''
    workdir = Path(workdir)
    workdir.mkdir(exist_ok=True, parents=True)
    cache_file = workdir / _cache_key(path_a, path_b)
    if cache_file.exists():
        with open(cache_file, 'rb') as f:
            return id_a, id_b, pickle.load(f)

    stem = cache_file.stem
    mat_file = workdir / f'mat_{stem}.txt'

    # -do is a boolean flag (no filename) - distances print inline in stdout
    result = subprocess.run(
        ['USalign', str(path_a), str(path_b), '-do', '-m', str(mat_file)],
        capture_output=True, text=True, timeout=120, check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f'USalign failed for {path_a} vs {path_b} (exit {result.returncode}):\n'
            f'--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}'
        )

    parsed = parse_usalign_output(result.stdout)
    parsed['rotation_matrix_file'] = str(mat_file) if mat_file.exists() else None

    with open(cache_file, 'wb') as f:
        pickle.dump(parsed, f)
    return id_a, id_b, parsed


def build_library(structures, workdir, max_workers=None):
    '''Run all pairwise USalign alignments in parallel and collect TM-scores.'''
    ids = list(structures)
    library, tm_scores = {}, {}

    with ProcessPoolExecutor(max_workers=max_workers) as pool:
        futures = {
            pool.submit(
                run_usalign_pair_cached, a, structures[a]['path'], b,
                structures[b]['path'], workdir): (a, b)
            for a, b in combinations(ids, 2)
        }
        for fut in as_completed(futures):
            id_a, id_b, parsed = fut.result()
            library[(id_a, id_b)] = parsed
            tm_scores[(id_a, id_b)] = 0.5 * (parsed['tm1'] + parsed['tm2'])
    return library, tm_scores



# guide tree
def build_guide_tree(ids, tm_scores):
    '''Build a hierarchical-clustering guide tree from pairwise TM-scores.'''
    n = len(ids)
    dist = np.zeros((n, n), dtype=np.float64)
    for i, a in enumerate(ids):
        for j, b in enumerate(ids):
            if i < j:
                key = (a, b) if (a, b) in tm_scores else (b, a)
                d = 1.0 - tm_scores[key]
                dist[i, j] = dist[j, i] = d
    return linkage(squareform(dist, checks=False), method='average')



# vectorized score matrix + Gotoh affine-gap NW
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


def build_msa(structures, library, tm_scores):
    '''Progressively build a multiple sequence alignment for a set of structures.'''
    ids = list(structures)
    n = len(ids)

    profiles = {i: [{ids[i]: r} for r in range(structures[ids[i]]['length'])] for i in range(n)}

    if n == 1:
        return profiles[0]
    if n == 2:
        return merge_profiles(profiles[0], profiles[1], library)

    guide_tree = build_guide_tree(ids, tm_scores)
    next_cluster_id = n
    for row in guide_tree:
        left, right = int(row[0]), int(row[1])
        profiles[next_cluster_id] = merge_profiles(profiles[left], profiles[right], library)
        del profiles[left], profiles[right]
        next_cluster_id += 1

    return profiles[next_cluster_id - 1]


def profile_to_aligned_sequences(root_profile, sequences):
    '''Convert a merged profile back into per-structure aligned sequence strings.'''
    aligned = {sid: [] for sid in sequences}
    for col in root_profile:
        for sid in sequences:
            aligned[sid].append(sequences[sid][col[sid]] if sid in col else '-')
    return {sid: ''.join(chars) for sid, chars in aligned.items()}


# Refinement
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





def _parse_args():
    '''Parse command-line arguments for the structure_msa CLI.'''
    parser = argparse.ArgumentParser(
        description='Structure-informed MSA using USalign + progressive merging'
    )
    parser.add_argument('--pdb-dir', type=Path, required=True, help='Directory of input PDB files')
    parser.add_argument(
        '--cache-dir', type=Path, default=Path('usalign_cache'),
        help='Cache directory for USalign outputs',
    )
    parser.add_argument(
        '--workers', type=int, default=None,
        help='Number of parallel processes (default: cpu count)',
    )
    parser.add_argument(
        '--refine', action='store_true',
        help='Run refinement passes after initial MSA'
        )
    parser.add_argument(
        '--pattern', type=str, default='*.pdb',
        help='Glob pattern for PDB files'
        )
    return parser.parse_args()


def _load_structures(pdb_dir, pattern):
    '''Parse sequences from all matching PDB files, skipping any with no CA atoms.'''
    pdb_files = sorted(pdb_dir.glob(pattern))
    if len(pdb_files) < 2:
        raise SystemExit(f'Need at least 2 PDB files in {pdb_dir}, found {len(pdb_files)}')

    print(f'Found {len(pdb_files)} structures:')
    structures, sequences = {}, {}
    for path in pdb_files:
        sid = path.stem
        seq = extract_sequence(path)
        if not seq:
            print(f'  WARNING: no CA atoms parsed from {path}, skipping')
            continue
        structures[sid] = {'path': path, 'length': len(seq)}
        sequences[sid] = seq
        print(f'  {sid}: {len(seq)} residues')

    if len(structures) < 2:
        raise SystemExit('Fewer than 2 structures parsed successfully, aborting.')
    return structures, sequences


def main():
    '''Run the structure_msa CLI end-to-end: align, build MSA, optionally refine, print.'''
    args = _parse_args()
    structures, sequences = _load_structures(args.pdb_dir, args.pattern)

    print('\nRunning pairwise USalign (cached, parallel)...')
    library, tm_scores = build_library(structures, args.cache_dir, args.workers)
    for (a, b), tm in sorted(tm_scores.items(), key=lambda x: -x[1]):
        print(f'  {a} vs {b}: TM-score = {tm:.3f}')

    print('\nBuilding progressive MSA...')
    root_profile = build_msa(structures, library, tm_scores)

    if args.refine:
        print('Refining...')
        root_profile = refine_msa(root_profile, library)

    aligned = profile_to_aligned_sequences(root_profile, sequences)

    print(f'\nFinal alignment ({len(root_profile)} columns):\n')
    max_id_len = max(len(sid) for sid in aligned)
    for sid, seq in aligned.items():
        print(f'{sid.ljust(max_id_len)}  {seq}')


if __name__ == '__main__':
    main()
