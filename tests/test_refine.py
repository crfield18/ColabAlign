'''The fast MSA refinement must make the same decisions as scoring whole profiles.'''

import numpy as np
import pytest

from colabalign.common import extract_sequence
from colabalign.consensus import (
    _column_luts, _resolve_pair_entry, build_cluster_msa, extract_leaf_profile, merge_profiles,
    refine_msa, remove_structure, structure_score,
)
from colabalign.usalign_runner import run_pairwise_alignments


def profile_score_sum(profile, library):
    '''Sum of the library weights the profile realises, over every pair of structures.'''
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


def _reference_refine(profile, library, max_passes=3):
    '''The previous implementation: re-score the whole profile for every candidate.'''
    struct_ids = sorted({sid for col in profile for sid in col})
    if len(struct_ids) < 3:
        return profile
    current, current_score = profile, profile_score_sum(profile, library)
    for _ in range(max_passes):
        improved = False
        for sid in struct_ids:
            candidate = merge_profiles(remove_structure(sid, current),
                                       extract_leaf_profile(sid, current), library)
            candidate_score = profile_score_sum(candidate, library)
            if candidate_score > current_score:
                current, current_score = candidate, candidate_score
                improved = True
        if not improved:
            break
    return current


def _cluster(structures, usalign, tmp_path):
    tm, store, _ = run_pairwise_alignments(structures, usalign, 2, tmp_path / 'cache')
    names = [p.stem for p in structures]
    sequences = {p.stem: extract_sequence(p) for p in structures}
    return names, sequences, store.sub_library(names), tm


@pytest.mark.parametrize('per_family', [3, 6])
def test_refine_matches_whole_profile_scoring(usalign, tmp_path, per_family):
    from conftest import make_structure_set
    structures = make_structure_set(tmp_path / 'm', per_family=per_family)
    names, sequences, library, tm = _cluster(structures, usalign, tmp_path)
    groups = [names] + [[n for n in names if n.startswith(f)] for f in ('hh', 'bab', 'hhh')]
    for members in groups:
        if len(members) < 3:
            continue
        seqs = {m: sequences[m] for m in members}
        start = build_cluster_msa(members, seqs, library, tm.loc[members, members])
        assert refine_msa(start, library) == _reference_refine(start, library)


def test_structure_score_is_the_profile_score_split_by_structure(usalign, tmp_path):
    from conftest import make_structure_set
    structures = make_structure_set(tmp_path / 'm', per_family=3)
    names, sequences, library, tm = _cluster(structures, usalign, tmp_path)
    profile = build_cluster_msa(names, sequences, library, tm)
    luts = _column_luts(profile)
    halves = sum(structure_score(s, luts, library, [o for o in names if o != s]) for s in names)
    assert np.isclose(halves / 2, profile_score_sum(profile, library))


def _reference_gotoh_matrices(score, gap_open=6.0, gap_extend=1.0):
    '''The previous cell-by-cell fill.'''
    from colabalign.consensus import NEG_INF
    n, m = score.shape
    match, gap_x, gap_y = (np.full((n + 1, m + 1), NEG_INF) for _ in range(3))
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
    return match, gap_x, gap_y


def _score_matrices():
    rng = np.random.default_rng(0)
    for n, m, zeros in [(1, 1, 0.0), (1, 7, 0.5), (9, 1, 0.5), (40, 55, 0.9), (70, 64, 0.97),
                        (150, 160, 0.98)]:
        score = rng.random((n, m)).astype(np.float32) * 3
        score[rng.random((n, m)) < zeros] = 0  # sparse scores give many exact ties
        yield score


def _reference_traceback(match, gap_x, gap_y):
    '''The previous Python traceback.'''
    i, j = match.shape[0] - 1, match.shape[1] - 1
    state = int(np.argmax([match[i, j], gap_x[i, j], gap_y[i, j]]))
    path = []
    while i > 0 or j > 0:
        if state == 0 and i > 0 and j > 0:
            path.append((i - 1, j - 1))
            i, j = i - 1, j - 1
            state = 0 if (i, j) == (0, 0) else int(
                np.argmax([match[i, j], gap_x[i, j], gap_y[i, j]]))
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


def test_compiled_fill_matches_cell_by_cell():
    from colabalign.consensus import NEG_INF, _fill_tables
    for score in _score_matrices():
        for new, old in zip(_fill_tables(score, 6.0, 1.0, float(NEG_INF)),
                            _reference_gotoh_matrices(score)):
            assert np.array_equal(new, old)


def test_gotoh_path_matches_python_traceback():
    from colabalign.consensus import gotoh_align
    for score in _score_matrices():
        assert gotoh_align(score) == _reference_traceback(*_reference_gotoh_matrices(score))


def test_compiles_without_cache_when_no_cache_folder_is_writable(monkeypatch):
    from colabalign import consensus
    calls = []

    def fake_njit(**options):
        calls.append(options)
        if options.get('cache'):
            raise RuntimeError('cannot cache function: no locator available')
        return lambda function: function

    monkeypatch.setattr(consensus, '_numba_njit', fake_njit)
    assert consensus.njit(len) is len
    assert calls == [{'cache': True, 'nogil': True}, {'nogil': True}]
