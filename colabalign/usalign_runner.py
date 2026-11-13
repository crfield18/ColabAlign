'''Runs all-vs-all USalign structural alignments in parallel.

Work is split into length-balanced batches ("tiles"). Chains are sorted by
residue count and cut into contiguous blocks of roughly equal total residues;
each tile is one block against itself or against a later block, so every
pair is aligned exactly once. USalign run time scales with L1 x L2, so tile
costs are estimated from residue counts (not file sizes, which are skewed by
hydrogens, alternate conformations and waters). Tiles run most expensive
first, a bounded number at a time, and workers return only TM-scores: the
full per-pair results go to the on-disk PairStore.
'''

import re
import subprocess
import tempfile
from collections import namedtuple
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from itertools import combinations, product
from math import ceil, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm.auto import tqdm as tqdm_auto

from .common import count_residues
from .pair_store import PairStore, file_hash

USALIGN_OPTS = ('-do', '-m', '-', '-het', '2')  # -het 2: also align MSE (see common.py)

# Fixed per-pair cost (process start-up, file parsing) in residue^2 units.
# Measured: ~2 ms per pair overhead vs ~2 us per residue^2 of alignment.
PAIR_OVERHEAD_COST = 1000
# Aim for at least this many tiles per core, so the pool can balance load...
TILES_PER_CORE = 8
# ...and at most roughly this many pairs per tile, so progress stays smooth.
TARGET_PAIRS_PER_TILE = 500
# Tiles whose longest chain is at most this long run as a single USalign
# process (-dir / -dir1 -dir2 list mode), which saves the per-pair start-up
# cost. That saving is only significant for small chains (~1.4x at <=110
# residues, ~1.05x at a median of 150). Set to 0 to always run per pair.
LIST_MODE_MAX_RESIDUES = 120
# Keep at most this many tiles queued per worker, so the parent never holds
# more than a handful of pending tasks or results.
IN_FLIGHT_PER_CORE = 2

Tile = namedtuple('Tile', 'a_start a_stop b_start b_stop n_pairs cost max_len')


DO_PATTERN = re.compile(
    r'^\s*CA\s+\w{3}\s+\S+\s+\d+\s+CA\s+\w{3}\s+\S+\s+\d+\s+([\d.]+)\s*$'
)
VALID_ALN_CHARS = set('ABCDEFGHIKLMNPQRSTVWYXZUO-')


def reverse_transformation_matrix(transform_mx_forward: np.ndarray) -> np.ndarray:
    '''Invert a forward USalign transform matrix to superpose in the opposite direction.'''
    assert transform_mx_forward.shape == (3, 4)
    translate_v_forward = transform_mx_forward[:, 0:1]
    rotation_mx_forward = transform_mx_forward[:, 1:4]
    rotation_mx_reverse = rotation_mx_forward.T
    translate_v_reverse = -np.dot(rotation_mx_reverse, translate_v_forward)
    return np.hstack((translate_v_reverse, rotation_mx_reverse))


def _extract_rotation_matrix(lines):
    '''Parse the 3x4 rotation/translation matrix block out of USalign stdout lines.'''
    start = None
    for idx, line in enumerate(lines):
        if 'rotation matrix' in line.lower():
            start = idx
            break
    if start is None:
        return None

    data_rows, row_idx = [], start + 2
    while len(data_rows) < 3 and row_idx < len(lines):
        parts = [p for p in lines[row_idx].split(' ') if p]
        if len(parts) >= 5:
            data_rows.append(parts[1:5])
        row_idx += 1
    if len(data_rows) < 3:
        return None
    return np.array(data_rows, dtype=np.float64)


def _compute_d0(length):
    '''TM-score-style d0, with the standard floor for short structures
    (avoids taking a fractional power of a negative number).'''
    return 0.5 if length < 15 else max(1.24 * (length - 15) ** (1 / 3) - 1.8, 0.5)


def _find_alignment_block(lines):
    '''Locate the (seq_a, match_line, seq_b) block within USalign stdout.'''
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped and set(stripped) <= VALID_ALN_CHARS and len(stripped) > 3:
            return i
    return None


def _build_residue_pairs(seq_a, match_line, seq_b):
    '''Walk the alignment block to build aligned residue-index pairs and weights.'''
    res_a_idx = res_b_idx = 0
    pairs_a, pairs_b, weights = [], [], []
    for ca, m, cb in zip(seq_a, match_line, seq_b):
        a_gap, b_gap = ca == '-', cb == '-'
        if not a_gap and not b_gap and m in (':', '.'):
            pairs_a.append(res_a_idx)
            pairs_b.append(res_b_idx)
            weights.append(1.0 if m == ':' else 0.5)
        if not a_gap:
            res_a_idx += 1
        if not b_gap:
            res_b_idx += 1
    return pairs_a, pairs_b, np.array(weights, dtype=np.float32)


def _weight_by_distance(weights, do_distances):
    '''Recompute weights from CA distances (TM-score-style falloff) when available.'''
    if len(do_distances) != len(weights) or len(weights) == 0:
        return weights
    d = np.array(do_distances, dtype=np.float32)
    d0 = _compute_d0(len(d))
    return 1.0 / (1.0 + (d / d0) ** 2)


def _parse_tm_scores(text):
    '''Extract the two directional TM-scores from USalign output text.'''
    tm_matches = re.findall(r'TM-score=\s*([\d.]+)', text)
    if len(tm_matches) < 2:
        raise ValueError(f'Could not parse TM-scores from USalign output:\n{text[:500]}')
    return float(tm_matches[0]), float(tm_matches[1])


def parse_usalign_stdout(stdout: bytes) -> dict:
    '''Parse default (full) USalign stdout: TM-scores, rotation matrix, and
    residue-level equivalences with -do-based confidence weights.'''
    text = stdout.decode('utf-8', errors='ignore')
    lines = text.splitlines()

    tm_1, tm_2 = _parse_tm_scores(text)

    aln_start = _find_alignment_block(lines)
    if aln_start is None or aln_start + 2 >= len(lines):
        raise ValueError(f'Could not locate alignment block in USalign output:\n{text[:500]}')

    seq_a, match_line, seq_b = lines[aln_start], lines[aln_start + 1], lines[aln_start + 2]
    do_distances = [float(m.group(1)) for l in lines if (m := DO_PATTERN.match(l))]

    pairs_a, pairs_b, weights = _build_residue_pairs(seq_a, match_line, seq_b)
    weights = _weight_by_distance(weights, do_distances)

    transform_mx_forward = _extract_rotation_matrix(lines)
    transform_mx_reverse = (
        reverse_transformation_matrix(transform_mx_forward)
        if transform_mx_forward is not None else None
    )

    return {
        'tm_forward': tm_1,
        'tm_reverse': tm_2,
        'transform_mx_forward': transform_mx_forward,
        'transform_mx_reverse': transform_mx_reverse,
        'res_a': np.array(pairs_a, dtype=np.int32),
        'res_b': np.array(pairs_b, dtype=np.int32),
        'weight': weights,
    }



# Batching (pure functions, unit-tested)

def canonical_order(stems, lengths):
    '''Indices sorted longest chain first (ties by name). This order defines
    each pair's orientation: the earlier structure is USalign's Structure_1.'''
    return sorted(range(len(stems)), key=lambda i: (-lengths[i], stems[i]))


def choose_block_count(n_models, cores, tiles_per_core=TILES_PER_CORE,
                       target_pairs_per_tile=TARGET_PAIRS_PER_TILE):
    '''Smallest block count k whose k(k+1)/2 tiles meet both batching targets.'''
    if n_models < 2:
        return 1
    n_pairs = n_models * (n_models - 1) // 2
    desired_tiles = max(tiles_per_core * max(int(cores), 1),
                        ceil(n_pairs / target_pairs_per_tile))
    k = ceil((sqrt(8 * desired_tiles + 1) - 1) / 2)
    return max(1, min(k, n_models))


def partition_blocks(lengths, n_blocks):
    '''
    Split lengths (already in canonical order) into at most n_blocks
    contiguous, non-empty blocks of roughly equal total residues.
    Returns a list of (start, stop) index ranges.
    '''
    n = len(lengths)
    if n == 0:
        return []
    n_blocks = max(1, min(int(n_blocks), n))
    weights = [max(int(length), 1) for length in lengths]
    total = sum(weights)
    bounds, start, acc = [], 0, 0
    for idx, weight in enumerate(weights):
        acc += weight
        blocks_left = n_blocks - len(bounds) - 1  # blocks still to open after this one
        items_left = n - idx - 1
        if blocks_left == 0 or items_left == 0:
            continue
        # cut at the next equal-residue boundary, or when every remaining
        # block needs exactly one of the remaining items
        if acc >= total * (len(bounds) + 1) / n_blocks or items_left == blocks_left:
            bounds.append((start, idx + 1))
            start = idx + 1
    bounds.append((start, n))
    return bounds


def build_tiles(blocks, lengths, overhead=PAIR_OVERHEAD_COST):
    '''
    One tile per block pair (a <= b), with estimated cost
    sum(L_i * L_j + overhead) over its pairs. Sorted most expensive first.
    '''
    sums = [sum(lengths[s:e]) for s, e in blocks]
    squares = [sum(x * x for x in lengths[s:e]) for s, e in blocks]
    tiles = []
    for a, (a_start, a_stop) in enumerate(blocks):
        for b in range(a, len(blocks)):
            b_start, b_stop = blocks[b]
            if a == b:
                size = a_stop - a_start
                n_pairs = size * (size - 1) // 2
                cost = (sums[a] ** 2 - squares[a]) / 2 + overhead * n_pairs
            else:
                n_pairs = (a_stop - a_start) * (b_stop - b_start)
                cost = sums[a] * sums[b] + overhead * n_pairs
            if n_pairs:
                # lengths are sorted descending, so block a's first entry is the longest
                tiles.append(Tile(a_start, a_stop, b_start, b_stop, n_pairs, cost,
                                  max(lengths[a_start], lengths[b_start])))
    tiles.sort(key=lambda t: t.cost, reverse=True)
    return tiles


def tile_pairs(tile):
    '''Canonical (i, j) index pairs, i < j, covered by a tile.'''
    if tile.a_start == tile.b_start:
        return combinations(range(tile.a_start, tile.a_stop), 2)
    return product(range(tile.a_start, tile.a_stop), range(tile.b_start, tile.b_stop))


# Workers

_WORKER = {}


def _init_worker(ctx):
    '''Receive the shared, read-only context once per worker process.'''
    _WORKER.clear()
    _WORKER.update(ctx)


def _run_usalign_pair(usalign_path, path_a, path_b):
    '''Align one pair. Never raises: returns (parsed_or_None, error_or_None).'''
    cmd = [str(usalign_path), str(path_a), str(path_b), *USALIGN_OPTS]
    try:
        result = subprocess.run(cmd, capture_output=True, check=False)
        if not result.stdout:
            stderr_text = result.stderr.decode(errors='ignore')
            return None, f'USalign produced no output (stderr: {stderr_text})'
        return parse_usalign_stdout(result.stdout), None
    except Exception as e:  # pylint: disable=broad-exception-caught
        # any failure (subprocess, parsing, etc.) is reported as a failed
        # pair rather than crashing the whole tile
        return None, f'{type(e).__name__}: {e}'


_NAME_PATTERNS = {
    k: re.compile(rf'^Name of Structure_{k}: (.+?)(?: \(to be superimposed.*)?$', re.M)
    for k in (1, 2)
}


def _structure_file_name(block, k):
    '''File name of Structure_k in one list-mode block ("dir/name.pdb:A" -> "name.pdb").'''
    match = _NAME_PATTERNS[k].search(block)
    if match is None:
        return None
    return Path(match.group(1).rsplit(':', 1)[0]).name


def split_list_mode_output(text):
    '''Split USalign -dir/-dir1/-dir2 stdout into one text block per aligned pair.'''
    return re.split(r'(?m)^(?=Name of Structure_1: )', text)[1:]


def _run_list_mode(usalign_path, paths, tile):
    '''
    Align a whole tile with a single USalign process. USalign holds one
    chain from list 1 and re-reads each list-2 chain from disk, so memory
    stays at two structures. Returns {(i, j): parsed} for the pairs it could
    parse; anything missing is retried per pair by the caller.
    '''
    members = list(range(tile.a_start, tile.a_stop)) + list(range(tile.b_start, tile.b_stop))
    index_by_name = {Path(paths[i]).name: i for i in members}
    models_dir = Path(paths[tile.a_start]).parent
    # list files are whitespace-separated names relative to one folder
    if (len(index_by_name) != len(set(members))
            or any(Path(paths[i]).parent != models_dir for i in members)
            or any(any(c.isspace() for c in name) for name in index_by_name)):
        return {}
    with tempfile.TemporaryDirectory() as tmp:
        list_a = Path(tmp, 'list_a')
        list_a.write_text(''.join(f'{Path(paths[i]).name}\n'
                                  for i in range(tile.a_start, tile.a_stop)),
                          encoding='utf-8')
        dir_arg = f'{models_dir.as_posix()}/'
        if tile.a_start == tile.b_start:
            cmd = [str(usalign_path), '-dir', dir_arg, str(list_a), *USALIGN_OPTS]
        else:
            list_b = Path(tmp, 'list_b')
            list_b.write_text(''.join(f'{Path(paths[j]).name}\n'
                                      for j in range(tile.b_start, tile.b_stop)),
                              encoding='utf-8')
            cmd = [str(usalign_path), '-dir1', dir_arg, str(list_a),
                   '-dir2', dir_arg, str(list_b), *USALIGN_OPTS]
        try:
            result = subprocess.run(cmd, capture_output=True, check=False)
        except OSError:
            return {}

    expected = set(tile_pairs(tile))
    parsed = {}
    for block in split_list_mode_output(result.stdout.decode('utf-8', errors='ignore')):
        i = index_by_name.get(_structure_file_name(block, 1))
        j = index_by_name.get(_structure_file_name(block, 2))
        if (i, j) not in expected:
            continue
        try:
            parsed[(i, j)] = parse_usalign_stdout(block.encode())
        except ValueError:
            continue  # retried per pair
    return parsed


def _run_tile(tile):
    '''
    Align every pair in a tile, reusing cached results. Saves full results
    to the PairStore and returns only [(i, j, tm_forward, tm_reverse, error)].
    '''
    paths, stems = _WORKER['paths'], _WORKER['stems']
    store, usalign_path = _WORKER['store'], _WORKER['usalign']

    results, todo = [], []
    for i, j in tile_pairs(tile):
        cached = store.load_canonical(stems[i], stems[j])
        if cached is not None:
            results.append((i, j, cached['tm_forward'], cached['tm_reverse'], None))
        else:
            todo.append((i, j))

    # List mode only for all-uncached tiles of small chains; it would redo cached pairs.
    if (len(todo) > 1 and len(todo) == tile.n_pairs
            and tile.max_len <= _WORKER['list_mode_max_residues']):
        for (i, j), entry in _run_list_mode(usalign_path, paths, tile).items():
            store.save(stems[i], stems[j], entry)
            results.append((i, j, entry['tm_forward'], entry['tm_reverse'], None))
        done = {(i, j) for i, j, *_ in results}
        todo = [pair for pair in todo if pair not in done]

    for i, j in todo:
        entry, error = _run_usalign_pair(usalign_path, paths[i], paths[j])
        if entry is None:
            results.append((i, j, None, None, error))
            continue
        store.save(stems[i], stems[j], entry)
        results.append((i, j, entry['tm_forward'], entry['tm_reverse'], None))
    return results


# Parent

# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def run_pairwise_alignments(model_list, usalign_path, cores, cache_dir,
                            list_mode_max_residues=LIST_MODE_MAX_RESIDUES,
                            tiles_per_core=TILES_PER_CORE,
                            target_pairs_per_tile=TARGET_PAIRS_PER_TILE):
    '''
    Run all-vs-all USalign in length-balanced tiles, with disk caching and
    per-pair error isolation.
    Returns: tmmatrix (pd.DataFrame, rows/cols in model_list order),
             pair_store (PairStore for transforms and residue equivalences),
             failed_pairs (list of (stem1, stem2, error)).
    '''
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True, parents=True)
    cores = max(int(cores), 1)

    names = [m.stem for m in model_list]
    lengths_by_name = [count_residues(m) for m in model_list]
    order = canonical_order(names, lengths_by_name)
    paths = [str(model_list[i]) for i in order]
    stems = [names[i] for i in order]
    lengths = [lengths_by_name[i] for i in order]

    store = PairStore(cache_dir,
                      hashes={stem: file_hash(path) for stem, path in zip(stems, paths)},
                      ranks={stem: rank for rank, stem in enumerate(stems)})

    blocks = partition_blocks(lengths, choose_block_count(
        len(stems), cores, tiles_per_core, target_pairs_per_tile))
    tiles = build_tiles(blocks, lengths)
    n_pairs = sum(t.n_pairs for t in tiles)
    print(f'Running US-align jobs for {len(stems)} chains: {n_pairs} pairs in '
          f'{len(tiles)} length-balanced batches.')

    position = {name: idx for idx, name in enumerate(names)}
    matrix_data = np.full((len(names), len(names)), np.nan)
    np.fill_diagonal(matrix_data, 1.0)
    failed_pairs = []

    total_bar = tqdm_auto(total=n_pairs, desc='Total alignments',
                          unit=' alignment', dynamic_ncols=True, leave=True)
    ctx = {'paths': paths, 'stems': stems, 'store': store, 'usalign': str(usalign_path),
           'list_mode_max_residues': list_mode_max_residues}

    with ProcessPoolExecutor(max_workers=cores, initializer=_init_worker,
                             initargs=(ctx,)) as executor:
        pending_tiles = iter(tiles)
        in_flight = {}

        def submit_next():
            tile = next(pending_tiles, None)
            if tile is not None:
                in_flight[executor.submit(_run_tile, tile)] = tile

        for _ in range(cores * IN_FLIGHT_PER_CORE):
            submit_next()

        while in_flight:
            done, _ = wait(in_flight, return_when=FIRST_COMPLETED)
            for future in done:
                tile = in_flight.pop(future)
                submit_next()
                try:
                    results = future.result()
                except Exception as e:  # pylint: disable=broad-exception-caught
                    # a lost tile (e.g. a crashed worker) only loses its own pairs
                    results = [(i, j, None, None, f'batch failed: {type(e).__name__}: {e}')
                               for i, j in tile_pairs(tile)]
                for i, j, tm_forward, tm_reverse, error in results:
                    if error is not None:
                        failed_pairs.append((stems[i], stems[j], error))
                        continue
                    a, b = position[stems[i]], position[stems[j]]
                    matrix_data[a, b] = matrix_data[b, a] = max(tm_forward, tm_reverse)
                total_bar.update(len(results))
    total_bar.close()

    if failed_pairs:
        print(f'\n{len(failed_pairs)} pair(s) failed and were excluded from the TM-score matrix:')
        for m1, m2, err in failed_pairs:
            print(f'  {m1} vs {m2}: {err}')

    tmmatrix = pd.DataFrame(matrix_data, index=names, columns=names)
    return tmmatrix, store, failed_pairs


def check_tm_matrix_completeness(tmmatrix):
    '''Returns list of (row, col) label pairs where the matrix has a NaN,
    excluding the diagonal. Call this before building the tree.'''
    missing = []
    n = tmmatrix.shape[0]
    labels = tmmatrix.index
    for i in range(n):
        for j in range(i + 1, n):
            if pd.isna(tmmatrix.iat[i, j]):
                missing.append((labels[i], labels[j]))
    return missing


