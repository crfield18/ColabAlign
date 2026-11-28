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


def _treecluster_groups(tree, threshold):
    '''(clusters with two or more members, singletons) from TreeCluster.py.'''
    out = subprocess.run(['TreeCluster.py', '-i', str(tree), '-t', format_threshold(threshold)],
                         capture_output=True, text=True, check=True).stdout.splitlines()[1:]
    labels = [line.split('\t')[1] for line in out]
    return len(set(labels) - {'-1'}), labels.count('-1')


@pytest.mark.skipif(shutil.which('TreeCluster.py') is None, reason='TreeCluster.py not on PATH')
@pytest.mark.parametrize('n,seed', [(12, 1), (40, 2), (75, 3)])
def test_cut_options_match_treecluster(tmp_path, n, seed):
    _, tree = _tree(tmp_path, n, seed)
    options = cut_options(tree, 0.0)
    groups = [o.n_groups for o in options]
    assert groups == sorted(groups, reverse=True) and len(set(groups)) == len(groups)
    assert groups[-1] == 1 and (options[-1].n_clusters, options[-1].n_singletons) == (1, 0)
    # the first cut offered has the most multi-member clusters
    assert options[0].n_clusters == max(o.n_clusters for o in options)
    rng = np.random.default_rng(seed)
    for option in rng.choice(options, size=min(12, len(options)), replace=False):
        assert option.threshold_min <= option.threshold < option.threshold_max
        assert _treecluster_groups(tree, option.threshold) == (option.n_clusters, option.n_singletons)


@pytest.mark.skipif(shutil.which('TreeCluster.py') is None, reason='TreeCluster.py not on PATH')
def test_cuts_below_the_peak_are_not_offered(tmp_path):
    _, tree = _tree(tmp_path, 40, 2)
    offered = cut_options(tree, 0.0)
    # just below the first offered cut, there are fewer multi-member clusters
    below = offered[0].threshold_min - 1e-9
    if below >= 0:
        assert _treecluster_groups(tree, below)[0] <= offered[0].n_clusters


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
    assert low_coverage_count(80, 0.2, 2) == 2           # only 2 structures left to add


def _options(*sizes):
    return [CutOption(k, 0, 0.1 * i, 0.1 * (i + 1), 0) for i, k in enumerate(sizes)]


def test_choose_cuts_closest_and_ranges():
    options = _options(30, 24, 20, 13, 9)
    assert [o.n_clusters for o in choose_cuts(options, [(20, 20)])] == [20]
    # 22 is 2 from 20 and 24: the smaller panel wins the tie
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(22, 22)])] == [20]
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(10, 25)])] == [24, 20, 13]
    # nothing in range -> closest to it
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(14, 18)])] == [13]
    assert [o.n_clusters for o in choose_cuts(_options(30, 24, 20, 13), [(100, 200)])] == [30]


def test_choose_cuts_include_additions():
    options = [CutOption(k, 50, 0.1 * i, 0.1 * (i + 1), low_coverage_count(k, 0.05, 50))
               for i, k in enumerate([60, 57, 50, 40])]
    chosen = choose_cuts(options, [(60, 60)])
    assert [(o.n_clusters, o.panel_size) for o in chosen] == [(57, 60)]


def test_threshold_is_short_and_inside_interval():
    assert CutOption(5, 0, 0.1534, 0.2011, 0).threshold == 0.18
    assert CutOption(5, 0, 0.15341, 0.15349, 0).threshold == 0.15345
    assert CutOption(1, 0, 0.734, float('inf'), 0).threshold == 0.74
    assert CutOption(9, 0, 0.0, 0.004, 0).threshold == 0.0


def test_add_low_coverage_is_farthest_first():
    names = ['a', 'b', 'c', 'd', 'e']
    tm = pd.DataFrame(np.array([
        [1.0, 0.9, 0.3, 0.3, 0.2],
        [0.9, 1.0, 0.3, 0.3, 0.2],
        [0.3, 0.3, 1.0, 0.95, 0.4],
        [0.3, 0.3, 0.95, 1.0, 0.4],
        [0.2, 0.2, 0.4, 0.4, 1.0],
    ]), index=names, columns=names)
    everyone = names
    # e is worst covered by {a}; after adding e, c/d tie at 0.4 -> c by name
    assert add_low_coverage(tm, ['a'], 2, everyone) == [('e', 0.2), ('c', 0.4)]
    # never more additions than candidates outside the panel
    assert [name for name, _ in add_low_coverage(tm, ['a', 'b'], 10, everyone)] == ['e', 'c', 'd']
    # only candidates are added, even when others are covered worse
    assert add_low_coverage(tm, ['a'], 5, ['c', 'd']) == [('c', 0.3), ('d', 0.95)]


def test_choose_cuts_same_size_keeps_fewest_singletons():
    # the singleton cap can give two cuts the same panel size: 10+2 and 12+0
    options = [CutOption(12, 0, 0.0, 0.1, 0), CutOption(10, 2, 0.1, 0.2, 2),
               CutOption(9, 1, 0.2, 0.3, 1)]
    assert [o.panel_size for o in options] == [12, 12, 10]
    assert [o.n_clusters for o in choose_cuts(options, [(12, 12)])] == [12]
    assert [o.n_clusters for o in choose_cuts(options, [(1, 20)])] == [12, 9]


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


def test_low_coverage_additions_off_by_default():
    from colabalign.cli import script_args
    assert script_args(['-i', 'x', '-o', 'y']).low_coverage_fraction == 0
    assert script_args(['-i', 'x', '-o', 'y', '-l', '0.05']).low_coverage_fraction == 0.05


def test_build_panel_takes_singletons_first_then_cluster_members(tmp_path):
    import csv
    from conftest import make_structure_set
    from colabalign.panel import build_panel
    files = make_structure_set(tmp_path / 'm', per_family=2)          # 6 structures
    names = sorted(p.stem for p in files)
    a, b, c, d, e, f = names
    tm = pd.DataFrame(0.5, index=names, columns=names)
    for x, y, score in [(a, b, 0.9), (a, c, 0.6), (d, e, 0.8), (f, a, 0.1), (f, d, 0.1),
                        (c, f, 0.1), (b, f, 0.1), (e, f, 0.1)]:
        tm.loc[x, y] = tm.loc[y, x] = score
    for name in names:
        tm.loc[name, name] = 1.0
    clusters = {1: [a, b, c], 2: [d, e], -1: [f]}
    reps = {1: {'record_name': a}, 2: {'record_name': d}}
    models = {p.stem: p for p in files}

    # 2 clusters with f = 0.5 -> 2 additions: the singleton, then the worst-covered member (c)
    row = build_panel('0.50', clusters, reps, tm, 0.5, models, tmp_path / 'out')
    with open(tmp_path / 'out' / 'threshold_0.50_panel.csv', newline='', encoding='utf-8') as h:
        panel = [(r['structure'], r['role'], r['cluster'], r['cluster_size']) for r in csv.DictReader(h)]
    assert panel == [(a, 'cluster_representative', '1', '3'), (d, 'cluster_representative', '2', '2'),
                     (f, 'low_coverage', '-1', '1'), (c, 'low_coverage', '1', '3')]
    assert (row['panel_size'], row['low_coverage_additions'], row['singletons'],
            row['singletons_in_panel']) == (4, 2, 1, 1)

    # with f = 0.2 -> 1 addition, and it is the singleton even though c is also poorly covered
    row = build_panel('0.20', clusters, reps, tm, 0.2, models, tmp_path / 'out2')
    with open(tmp_path / 'out2' / 'threshold_0.20_panel.csv', newline='', encoding='utf-8') as h:
        assert [r['structure'] for r in csv.DictReader(h)] == [a, d, f]
