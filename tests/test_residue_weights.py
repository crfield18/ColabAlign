"""Regressions for shorter-full-chain distance scaling in the structure MSA."""

import numpy as np
import pytest

from colabalign.usalign_runner import (
    _compute_d0,
    _parse_structure_lengths,
    _weight_by_distance,
    parse_usalign_stdout,
)


def _stdout(length_a=200, length_b=200, aligned=50, distance=3.0):
    """Full sequences include unaligned termini, just like USalign stdout."""
    sequence_a = "A" * length_a + "-" * (length_b - aligned)
    sequence_b = "-" * (length_a - aligned) + "G" * length_b
    matches = " " * (length_a - aligned) + ":" * aligned + " " * (length_b - aligned)
    distances = "\n".join(
        f" CA ALA A {length_a - aligned + index + 1} CA GLY B {index + 1} {distance:.3f}"
        for index in range(aligned)
    )
    return (
        f"Name of Structure_1: a.pdb:A\n"
        f"Name of Structure_2: b.pdb:B\n"
        f"Length of Structure_1: {length_a} residues\n"
        f"Length of Structure_2: {length_b} residues\n"
        f"Aligned length= {aligned}, RMSD= 3.00\n"
        f"TM-score= 0.40000 (normalized by length of Structure_1)\n"
        f"TM-score= 0.50000 (normalized by length of Structure_2)\n"
        f"{sequence_a}\n{matches}\n{sequence_b}\n"
        f"{distances}\n"
    ).encode()


@pytest.mark.parametrize("length", [1, 8, 14, 15, 20, 21])
def test_short_chain_d0_floor(length):
    assert _compute_d0(length) == 0.5


@pytest.mark.parametrize(
    "length, expected",
    [(22, 0.572034), (50, 2.256123), (100, 3.652069), (200, 5.265544)],
)
def test_full_chain_d0_formula(length, expected):
    assert _compute_d0(length) == pytest.approx(expected, abs=1e-5)


def test_weights_are_not_divided_by_chain_length():
    weights = _weight_by_distance(np.ones(3), [0.0, 3.0, 6.0], 200)
    assert weights[0] == 1.0
    assert weights[1] == pytest.approx(0.75494, abs=1e-5)
    assert np.all((weights > 0) & (weights <= 1))
    assert weights[0] > weights[1] > weights[2]


def test_weights_do_not_depend_on_aligned_pair_count():
    two = _weight_by_distance(np.ones(2), [3.0, 6.0], 200)
    three = _weight_by_distance(np.ones(3), [3.0, 6.0, 2.0], 200)
    np.testing.assert_array_equal(two, three[:2])


@pytest.mark.parametrize("length_a,length_b", [(200, 150), (150, 200)])
def test_parser_uses_shorter_full_chain_in_either_orientation(length_a, length_b):
    entry = parse_usalign_stdout(_stdout(length_a, length_b, aligned=25))
    assert entry["weight_length"] == 150
    np.testing.assert_allclose(entry["weight"], 0.69802, atol=1e-5)
    assert entry["tm_forward"] == 0.4
    assert entry["tm_reverse"] == 0.5
    np.testing.assert_array_equal(entry["res_a"], np.arange(length_a - 25, length_a))
    np.testing.assert_array_equal(entry["res_b"], np.arange(25))


def test_partial_alignment_uses_full_length_not_aligned_region():
    entry = parse_usalign_stdout(_stdout())
    assert len(entry["res_a"]) == len(entry["res_b"]) == 50
    assert entry["weight_length"] == 200
    np.testing.assert_allclose(entry["weight"], 0.75494, atol=1e-5)
    old_aligned_length_weight = 1.0 / (1.0 + (3.0 / _compute_d0(50)) ** 2)
    assert not np.allclose(entry["weight"], old_aligned_length_weight)


def test_parser_weights_are_invariant_to_coverage_for_fixed_chain_lengths():
    partial = parse_usalign_stdout(_stdout(aligned=25))
    full = parse_usalign_stdout(_stdout(aligned=200))
    np.testing.assert_array_equal(partial["weight"], full["weight"][:25])


def test_parser_uses_standard_short_chain_floor():
    entry = parse_usalign_stdout(_stdout(length_a=10, length_b=8, aligned=4))
    assert entry["weight_length"] == 8
    np.testing.assert_allclose(entry["weight"], 1.0 / 37.0)


@pytest.mark.parametrize("distances", [[], [3.0]])
def test_unavailable_or_mismatched_distances_preserve_match_symbol_weights(distances):
    original = np.array([1.0, 0.5], dtype=np.float32)
    assert _weight_by_distance(original, distances, 200) is original


def test_no_aligned_pairs_keeps_empty_weights():
    original = np.array([], dtype=np.float32)
    assert _weight_by_distance(original, [], 200) is original


@pytest.mark.parametrize("length", [0, -1])
def test_invalid_normalizing_length_is_rejected(length):
    with pytest.raises(ValueError, match="must be positive"):
        _weight_by_distance(np.ones(1), [3.0], length)


@pytest.mark.parametrize(
    "text",
    [
        "",
        "Length of Structure_1: 200 residues\n",
        "Length of Structure_1: 0 residues\nLength of Structure_2: 50 residues\n",
    ],
)
def test_missing_or_invalid_full_lengths_do_not_revert_to_aligned_pair_count(text):
    with pytest.raises(ValueError, match="positive full chain lengths"):
        _parse_structure_lengths(text)


def test_parser_requires_both_reported_full_lengths():
    stdout = _stdout().replace(b"Length of Structure_2: 200 residues\n", b"")
    with pytest.raises(ValueError, match="positive full chain lengths"):
        parse_usalign_stdout(stdout)


@pytest.mark.parametrize("missing_distances", [1, 50])
def test_parser_incomplete_distance_table_uses_symbol_fallback(missing_distances):
    lines = _stdout().decode().splitlines()
    distance_indices = [index for index, line in enumerate(lines) if line.startswith(" CA")]
    remove = set(distance_indices[:missing_distances])
    stdout = "\n".join(line for index, line in enumerate(lines) if index not in remove).encode()
    entry = parse_usalign_stdout(stdout)
    assert entry["weight_length"] == 200
    np.testing.assert_array_equal(entry["weight"], np.ones(50))
    assert (entry["tm_forward"], entry["tm_reverse"]) == (0.4, 0.5)
