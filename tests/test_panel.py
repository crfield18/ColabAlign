'''Panel-size targeting, threshold formatting and low-coverage additions.'''

import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from colabalign.dendrogram import grow_tree
from colabalign.panel import (
    CutOption, add_low_coverage, choose_cuts, cut_options, format_threshold,
    low_coverage_count, parse_panel_size,
)


def _tm_matrix(n, seed):
    '''A symmetric TM-like matrix with loose groups, values in (0, 1], diagonal 1.'''
    rng = np.random.default_rng(seed)
    centres = rng.normal(scale=3, size=(max(2, n // 6), 4))
    points = centres[rng.integers(len(centres), size=n)] + rng.normal(scale=0.8, size=(n, 4))
    dist = np.linalg.norm(points[:, None] - points[None], axis=2)
    tm = np.exp(-dist / 4)
    np.fill_diagonal(tm, 1.0)
    names = [f's{i:03d}' for i in range(n)]
    return pd.DataFrame(tm, index=names, columns=names)


def _tree(tmp_path, n, seed):
    tm = _tm_matrix(n, seed)
    grow_tree(tm, [Path(f'{name}.pdb') for name in tm.index], tmp_path)
    return tm, tmp_path / 'colabalign.tree'


def _treecluster_count(tree, threshold):
    out = subprocess.run(['TreeCluster.py', '-i', str(tree), '-t', format_threshold(threshold)],
                         capture_output=True, text=True, check=True).stdout.splitlines()[1:]
    labels = [line.split('\t')[1] for line in out]
    return labels.count('-1') + len(set(labels) - {'-1'})


@pytest.mark.skipif(shutil.which('TreeCluster.py') is None, reason='TreeCluster.py not on PATH')
@pytest.mark.parametrize('n,seed', [(12, 1), (40, 2), (75, 3)])
def test_cut_options_match_treecluster(tmp_path, n, seed):
    _, tree = _tree(tmp_path, n, seed)
    options = cut_options(tree, 0.0)
    counts = [o.n_clusters for o in options]
    assert counts == sorted(counts, reverse=True) and len(set(counts)) == len(counts)
    assert counts[0] <= n and counts[-1] == 1
    rng = np.random.default_rng(seed)
    for option in rng.choice(options, size=min(12, len(options)), replace=False):
        assert option.threshold_min <= option.threshold < option.threshold_max
        assert _treecluster_count(tree, option.threshold) == option.n_clusters


def test_parse_panel_size():
    assert parse_panel_size('50') == (50, 50)
    assert parse_panel_size(' 40-60 ') == (40, 60)
    for bad in ('0', '60-40', 'abc', '1-2-3', '-5', '4.5'):
        with pytest.raises(ValueError):
            parse_panel_size(bad)


def test_format_threshold():
    assert format_threshold(0.2) == '0.20'
    assert format_threshold(0.155) == '0.155'
    assert format_threshold(0.1234) == '0.1234'
    assert format_threshold(1) == '1.00'


def test_low_coverage_count():
    # additions make up the fraction of the final panel, rounded half up
    assert low_coverage_count(19, 0.05, 100) == 1        # 1 of 20
    assert low_coverage_count(57, 0.05, 100) == 3        # 3 of 60
    assert low_coverage_count(9, 0.05, 100) == 0         # 0.47 rounds down
    assert low_coverage_count(10, 0.05, 100) == 1        # 0.53 rounds up
    assert low_coverage_count(40, 0.0, 100) == 0
    assert low_coverage_count(98, 0.2, 100) == 2         # cannot exceed the pool


def _options(*sizes):
    return [CutOption(k, 0.1 * i, 0.1 * (i + 1), 0) for i, k in enumerate(sizes)]


def test_choose_cuts_closest_and_ranges():
    options = _options(30, 24, 20, 13, 9)
    assert [o.n_clusters for o in choose_cuts(options, [(20, 20)])] == [20]
    # 22 is 2 from 20 and 24: the smaller panel wins the tie
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(22, 22)])] == [20]
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(10, 25)])] == [13, 20, 24]
    # nothing in range -> closest to it
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(14, 18)])] == [13]
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(100, 200)])] == [30]


def test_choose_cuts_include_additions():
    options = [CutOption(k, 0.1 * i, 0.1 * (i + 1), low_coverage_count(k, 0.05, 500))
               for i, k in enumerate([60, 57, 50, 40])]
    chosen = choose_cuts(options, [(60, 60)])
    assert [(o.n_clusters, o.panel_size) for o in chosen] == [(57, 60)]


def test_threshold_is_short_and_inside_interval():
    assert CutOption(5, 0.1534, 0.2011, 0).threshold == 0.18
    assert CutOption(5, 0.15341, 0.15349, 0).threshold == 0.15345
    assert CutOption(1, 0.734, float('inf'), 0).threshold == 0.74
    assert CutOption(9, 0.0, 0.004, 0).threshold == 0.0


def test_add_low_coverage_is_farthest_first():
    names = ['a', 'b', 'c', 'd', 'e']
    tm = pd.DataFrame(np.array([
        [1.0, 0.9, 0.3, 0.3, 0.2],
        [0.9, 1.0, 0.3, 0.3, 0.2],
        [0.3, 0.3, 1.0, 0.95, 0.4],
        [0.3, 0.3, 0.95, 1.0, 0.4],
        [0.2, 0.2, 0.4, 0.4, 1.0],
    ]), index=names, columns=names)
    # e is worst covered by {a}; after adding e, c/d tie at 0.4 -> c by name
    assert add_low_coverage(tm, ['a'], 2) == [('e', 0.2), ('c', 0.4)]
    # never more additions than structures outside the panel
    assert [name for name, _ in add_low_coverage(tm, ['a', 'b'], 10)] == ['e', 'c', 'd']


def test_choose_cuts_same_size_keeps_most_clusters():
    # 12 structures with 20 % additions: the pool-size cap gives 12 = 12+0 = 11+1 = 10+2
    options = [CutOption(k, 0.1 * i, 0.1 * (i + 1), low_coverage_count(k, 0.2, 12))
               for i, k in enumerate([12, 11, 10, 9])]
    assert [o.panel_size for o in options] == [12, 12, 12, 11]
    assert [o.n_clusters for o in choose_cuts(options, [(12, 12)])] == [12]
    assert [o.n_clusters for o in choose_cuts(options, [(1, 20)])] == [9, 12]


def test_tree_is_reused_only_for_identical_inputs(tmp_path, capsys):
    tm = _tm_matrix(15, 4)
    models = [Path(f'{name}.pdb') for name in tm.index]
    grow_tree(tm, models, tmp_path)
    first = (tmp_path / 'colabalign.tree').read_text()
    grow_tree(tm, models, tmp_path)
    assert 'Reusing the existing dendrogram' in capsys.readouterr().out
    assert (tmp_path / 'colabalign.tree').read_text() == first
    changed = tm.copy()
    changed.iloc[0, 1] = changed.iloc[1, 0] = 0.01
    grow_tree(changed, models, tmp_path)
    assert 'Generating structural dendrogram' in capsys.readouterr().out
    assert (tmp_path / 'colabalign.tree').read_text() != first
