'''Pure tests for the length-balanced tile batching.'''

from itertools import combinations

import pytest

from colabalign.usalign_runner import (
    Tile, build_tiles, canonical_order, choose_block_count, partition_blocks, tile_pairs,
)


def _lengths(n, seed=0):
    import random
    rng = random.Random(seed)
    return sorted((rng.randint(30, 900) for _ in range(n)), reverse=True)


@pytest.mark.parametrize('n,k', [(1, 1), (2, 1), (5, 3), (10, 10), (37, 6), (200, 17)])
def test_partition_blocks_contiguous_and_non_empty(n, k):
    lengths = _lengths(n)
    blocks = partition_blocks(lengths, k)
    assert blocks[0][0] == 0 and blocks[-1][1] == n
    assert all(start < stop for start, stop in blocks)
    assert all(a[1] == b[0] for a, b in zip(blocks, blocks[1:]))
    assert len(blocks) == min(k, n)


def test_partition_blocks_balances_residues():
    lengths = _lengths(400, seed=3)
    blocks = partition_blocks(lengths, 10)
    sums = [sum(lengths[s:e]) for s, e in blocks]
    # each block within one maximum chain length of the ideal share
    ideal = sum(lengths) / len(blocks)
    assert max(abs(s - ideal) for s in sums) <= max(lengths)
    # long chains end up in small blocks
    sizes = [e - s for s, e in blocks]
    assert sizes[0] < sizes[-1]


def test_partition_blocks_empty_and_oversized():
    assert partition_blocks([], 4) == []
    assert partition_blocks([50, 40], 10) == [(0, 1), (1, 2)]


def test_choose_block_count():
    assert choose_block_count(1, 8) == 1
    assert choose_block_count(2, 1) <= 2
    k = choose_block_count(1000, 4)
    n_tiles = k * (k + 1) // 2
    assert n_tiles >= 4 * 8                       # enough tiles per core
    assert 1000 * 999 / 2 / n_tiles <= 500         # tiles are not too large
    assert choose_block_count(5, 64) == 5          # never more blocks than models


@pytest.mark.parametrize('n,k', [(2, 1), (7, 3), (30, 6), (101, 13)])
def test_tiles_cover_every_pair_exactly_once(n, k):
    lengths = _lengths(n, seed=n)
    tiles = build_tiles(partition_blocks(lengths, k), lengths)
    pairs = [p for t in tiles for p in tile_pairs(t)]
    assert len(pairs) == len(set(pairs)) == n * (n - 1) // 2
    assert set(pairs) == set(combinations(range(n), 2))
    assert all(i < j for i, j in pairs)
    assert sum(t.n_pairs for t in tiles) == len(pairs)
    assert [t.cost for t in tiles] == sorted((t.cost for t in tiles), reverse=True)


def test_tile_cost_matches_pairwise_sum():
    lengths = _lengths(23, seed=9)
    for tile in build_tiles(partition_blocks(lengths, 4), lengths, overhead=7):
        expected = sum(lengths[i] * lengths[j] + 7 for i, j in tile_pairs(tile))
        assert tile.cost == pytest.approx(expected)
        assert tile.max_len == max(max(lengths[i], lengths[j]) for i, j in tile_pairs(tile))


def test_tiles_are_more_even_than_equal_count_chunks():
    lengths = _lengths(300, seed=1)
    tiles = build_tiles(partition_blocks(lengths, 12), lengths)
    off_diag = [t.cost for t in tiles if t.a_start != t.b_start]
    # equal-residue blocks make every off-diagonal tile cost about the same
    assert max(off_diag) / min(off_diag) < 2.0


def test_canonical_order_longest_first_ties_by_name():
    stems = ['b', 'a', 'c', 'd']
    lengths = [100, 100, 300, 50]
    assert [stems[i] for i in canonical_order(stems, lengths)] == ['c', 'a', 'b', 'd']


def test_tile_pairs_diagonal_and_rectangle():
    assert list(tile_pairs(Tile(0, 3, 0, 3, 3, 0, 0))) == [(0, 1), (0, 2), (1, 2)]
    assert list(tile_pairs(Tile(0, 2, 2, 4, 4, 0, 0))) == [(0, 2), (0, 3), (1, 2), (1, 3)]
