'''Runs pairwise USalign structural alignments in parallel, with caching and error isolation.'''

import hashlib
import pickle
import re
import subprocess
from itertools import combinations
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm.auto import tqdm as tqdm_auto

DO_PATTERN = re.compile(
    r'^\s*CA\s+\w{3}\s+\S+\s+\d+\s+CA\s+\w{3}\s+\S+\s+\d+\s+([\d.]+)\s*$'
)
VALID_ALN_CHARS = set('ABCDEFGHIKLMNPQRSTVWYXZUO-')


def _chunks(lst, n):
    '''Yield n chunks from the list.'''
    avg = len(lst) / n
    last = 0
    while last < len(lst):
        yield lst[int(last):int(last + avg)]
        last += avg


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


def _file_hash(path):
    '''Return a short sha256 hex digest of a file's contents, for cache keys.'''
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def _cache_key(path_a, path_b):
    '''Build a cache filename from the hashes of two input structure files.'''
    return f'{_file_hash(path_a)}_{_file_hash(path_b)}.pkl'


def _load_cached_pair(cache_file):
    '''Load a cached alignment result, returning None if missing or corrupt.'''
    if not cache_file.exists():
        return None
    try:
        with open(cache_file, 'rb') as f:
            return pickle.load(f)
    except Exception:  # pylint: disable=broad-exception-caught
        # any corrupt/unreadable cache entry should just trigger a recompute
        return None


def _run_usalign_pair(usalign_path, models_path, cache_dir, model1, model2):
    '''
    Never raises - returns (model1, model2, parsed_dict_or_None, error_or_None)
    so that one bad pair can't take down the rest of its chunk's results.
    '''
    path_a = models_path.joinpath(model1.name)
    path_b = models_path.joinpath(model2.name)

    cache_file = cache_dir.joinpath(_cache_key(path_a, path_b))
    cached = _load_cached_pair(cache_file)
    if cached is not None:
        return model1, model2, cached, None

    cmd = [
        str(usalign_path), str(path_a), str(path_b),
        '-do', '-m', '-', '-het', '2',  # -het 2: also align MSE (see common.py)
    ]
    try:
        with subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
            stdout, stderr = process.communicate()
        if not stdout:
            return model1, model2, None, f'USalign produced no output (stderr: {stderr.decode(errors='ignore')})'
        parsed = parse_usalign_stdout(stdout)
    except Exception as e:  # pylint: disable=broad-exception-caught
        # any failure (subprocess, parsing, etc.) is reported as a failed
        # pair rather than crashing the whole run
        return model1, model2, None, f'{type(e).__name__}: {e}'

    cache_dir.mkdir(exist_ok=True, parents=True)
    with open(cache_file, 'wb') as f:
        pickle.dump(parsed, f)

    return model1, model2, parsed, None


def _usalign_job(usalign_path, models_path, cache_dir, pairs):
    '''Run USalign for every pair in a chunk, sequentially, within one worker process.'''
    return [_run_usalign_pair(usalign_path, models_path, cache_dir, m1, m2) for m1, m2 in pairs]


def _dispatch_jobs(executor, usalign_path, models_path, cache_dir, jobs, total_bar):
    '''Submit all chunked jobs to the executor, wiring up progress-bar callbacks.'''
    futures = []
    for job in jobs:
        fut = executor.submit(_usalign_job, usalign_path, models_path, cache_dir, job)

        def _callback(f):
            try:
                n = len(f.result())
            except Exception:  # pylint: disable=broad-exception-caught
                # a failed chunk shouldn't stall the progress bar
                n = 0
            total_bar.update(n)

        fut.add_done_callback(_callback)
        futures.append(fut)
    return futures


def _collect_results(futures):
    '''Gather per-pair results from completed futures, tracking chunk-level failures.'''
    usalign_results_map = defaultdict(dict)
    failed_pairs = []

    for future in as_completed(futures):
        try:
            results = future.result()
        except Exception as e:  # pylint: disable=broad-exception-caught
            # one lost chunk shouldn't abort collection of the rest
            print(f'Chunk-level error (entire chunk lost): {e}')
            continue

        for model1, model2, parsed, error in results:
            if error is not None or parsed is None:
                failed_pairs.append((model1.stem, model2.stem, error))
                print(f'Skipping pair {model1.stem} vs {model2.stem}: {error}')
                continue
            usalign_results_map[model1.stem][model2.stem] = parsed

    return usalign_results_map, failed_pairs


def run_pairwise_alignments(model_list, models_path, usalign_path, cores, cache_dir):
    '''
    Runs all-pairs USalign in parallel, with disk caching and per-pair error
    isolation.
    Returns: tmmatrix (pd.DataFrame), library (dict for MSA scoring),
             transforms (dict for cluster superposition), failed_pairs (list)
    '''
    cache_dir.mkdir(exist_ok=True, parents=True)
    all_combos = list(combinations(model_list, 2))
    print(f'Running US-align jobs for {len(model_list)} chains.')

    chunk_count = int(cores) if len(all_combos) <= int(cores * 5000) else int(cores * 5000)
    jobs = list(_chunks(all_combos, chunk_count))

    total_bar = tqdm_auto(total=len(all_combos), desc='Total alignments',
                           unit=' alignment', dynamic_ncols=True, leave=True)

    with ProcessPoolExecutor(max_workers=cores) as executor:
        futures = _dispatch_jobs(executor, usalign_path, models_path, cache_dir, jobs, total_bar)
        usalign_results_map, failed_pairs = _collect_results(futures)

    if failed_pairs:
        print(f'\n{len(failed_pairs)} pair(s) failed and were excluded from the TM-score matrix:')
        for m1, m2, err in failed_pairs:
            print(f'  {m1} vs {m2}: {err}')

    tmmatrix = _build_tm_matrix(model_list, usalign_results_map)
    library = _build_library(usalign_results_map)
    transforms = _build_transforms(usalign_results_map)

    return tmmatrix, library, transforms, failed_pairs


def _build_tm_matrix(model_list, usalign_results_map):
    '''Assemble the symmetric TM-score matrix (max of forward/reverse) as a DataFrame.'''
    matrix_data = np.full((len(model_list), len(model_list)), np.nan)
    model_index = {model.stem: idx for idx, model in enumerate(model_list)}

    for m1, subsubmap in usalign_results_map.items():
        for m2, entry in subsubmap.items():
            max_tm = max(entry['tm_forward'], entry['tm_reverse'])
            i, j = model_index[m1], model_index[m2]
            matrix_data[i, j] = max_tm
            matrix_data[j, i] = max_tm

    np.fill_diagonal(matrix_data, 1.0000)

    return pd.DataFrame(
        matrix_data,
        index=[model.stem for model in model_list],
        columns=[model.stem for model in model_list]
    )


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


def _build_library(usalign_results_map):
    '''Build the residue-pair alignment library used for MSA profile scoring.'''
    library = {}
    for m1, subsubmap in usalign_results_map.items():
        for m2, entry in subsubmap.items():
            library[(m1, m2)] = {
                'res_a': entry['res_a'],
                'res_b': entry['res_b'],
                'weight': entry['weight'],
            }
    return library


def _build_transforms(usalign_results_map):
    '''Build the forward/reverse transform-matrix lookup used for cluster superposition.'''
    transforms = {}
    for m1, subsubmap in usalign_results_map.items():
        for m2, entry in subsubmap.items():
            transforms[(m1, m2)] = {
                'forward': entry['transform_mx_forward'],
                'reverse': entry['transform_mx_reverse'],
            }
    return transforms


def get_transform_matrix(transforms, source_name, target_name):
    '''
    Returns the matrix that superposes source_name onto target_name's frame,
    or None if the pair isn't in the library.
    '''
    key = (source_name, target_name)
    if key in transforms:
        return transforms[key]['forward']
    key = (target_name, source_name)
    if key in transforms:
        return transforms[key]['reverse']
    return None
