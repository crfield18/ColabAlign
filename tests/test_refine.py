'''The fast MSA refinement must make the same decisions as scoring whole profiles.'''

import numpy as np
import pytest

from colabalign.common import extract_sequence
from colabalign.consensus import (
    _column_luts, build_cluster_msa, extract_leaf_profile, merge_profiles,
    profile_score_sum, refine_msa, remove_structure, structure_score,
)
from colabalign.usalign_runner import run_pairwise_alignments


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


def test_gotoh_matches_cell_by_cell(monkeypatch):
    import colabalign.consensus as consensus
    rng = np.random.default_rng(0)
    captured = {}
    real_full = np.full

    def recording_full(shape, value):
        array = real_full(shape, value)
        captured.setdefault('arrays', []).append(array)
        return array

    for n, m, zeros in [(1, 1, 0.0), (1, 7, 0.5), (9, 1, 0.5), (40, 55, 0.9), (70, 64, 0.97)]:
        score = rng.random((n, m)).astype(np.float32) * 3
        score[rng.random((n, m)) < zeros] = 0  # sparse scores give many exact ties
        captured.clear()
        monkeypatch.setattr(consensus.np, 'full', recording_full)
        consensus.gotoh_align(score)
        monkeypatch.setattr(consensus.np, 'full', real_full)
        for new, old in zip(captured['arrays'], _reference_gotoh_matrices(score)):
            assert np.array_equal(new, old)
