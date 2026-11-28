'''Integration tests against the real external tools (skipped when not installed).'''

import stat

import numpy as np
import pandas as pd
import pytest

from colabalign.common import count_residues
from colabalign.usalign_runner import run_pairwise_alignments

REAL = {'list_mode_max_residues': 10_000}     # always use USalign's -dir list mode
PAIRWISE = {'list_mode_max_residues': 0}      # always run one USalign call per pair


def _run(structures, usalign, cache, cores=2, **kw):
    return run_pairwise_alignments(structures, usalign, cores, cache, **kw)


def _fake_usalign(tmp_path, real, body):
    script = tmp_path / 'fake_usalign'
    script.write_text(f'#!/bin/sh\n{body}\nexec "{real}" "$@"\n')
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def test_matrix_shape_and_values(structures, usalign, tmp_path):
    tm, store, failed = _run(structures, usalign, tmp_path / 'c')
    names = [p.stem for p in structures]
    assert not failed
    assert list(tm.index) == list(tm.columns) == names
    assert np.allclose(np.diag(tm.values), 1.0)
    assert np.allclose(tm.values, tm.values.T)
    assert not tm.isna().any().any()
    assert ((tm.values > 0) & (tm.values <= 1)).all()
    # same-family pairs should score higher on average than cross-family pairs
    fam = [n.rstrip('0123456789') for n in names]
    same = [tm.iloc[i, j] for i in range(len(names)) for j in range(i) if fam[i] == fam[j]]
    diff = [tm.iloc[i, j] for i in range(len(names)) for j in range(i) if fam[i] != fam[j]]
    assert np.mean(same) > np.mean(diff)
    n = len(names)
    assert len(list(store.cache_dir.glob('*.pkl'))) == n * (n - 1) // 2


def test_list_mode_matches_per_pair(structures, usalign, tmp_path):
    a, sa, _ = _run(structures, usalign, tmp_path / 'list', **REAL)
    b, sb, _ = _run(structures, usalign, tmp_path / 'pair', **PAIRWISE)
    pd.testing.assert_frame_equal(a, b)
    names = [p.stem for p in structures]
    for x, y in [(names[0], names[-1]), (names[-1], names[0]), (names[2], names[5])]:
        assert np.allclose(sa.transform(x, y), sb.transform(x, y))
    la, lb = sa.sub_library(names), sb.sub_library(names)
    assert la.keys() == lb.keys()
    by_name = {path.stem: path for path in structures}
    for key, entry in la.items():
        for field, value in entry.items():
            assert np.array_equal(value, lb[key][field])
        expected_length = min(count_residues(by_name[name]) for name in key)
        for store in (sa, sb):
            cached, _ = store.load(*key)
            assert cached['weight_length'] == expected_length


@pytest.mark.parametrize('cores,tiles_per_core,target', [(1, 1, 10_000), (3, 8, 5)])
def test_batching_does_not_change_results(structures, usalign, tmp_path, cores, tiles_per_core,
                                          target):
    ref, _, _ = _run(structures, usalign, tmp_path / 'ref', cores=1, tiles_per_core=1,
                     target_pairs_per_tile=10_000)
    got, _, _ = _run(structures, usalign, tmp_path / f'b{cores}', cores=cores,
                     tiles_per_core=tiles_per_core, target_pairs_per_tile=target)
    pd.testing.assert_frame_equal(ref, got)


def test_cache_is_reused(structures, usalign, tmp_path):
    first, _, _ = _run(structures, usalign, tmp_path / 'c')
    # a non-existent USalign proves nothing is recomputed
    second, _, failed = _run(structures, tmp_path / 'no_such_usalign', tmp_path / 'c')
    assert not failed
    pd.testing.assert_frame_equal(first, second)


def test_cache_survives_input_order_and_renames(structures, usalign, tmp_path):
    first, _, _ = _run(structures, usalign, tmp_path / 'c')
    renamed = []
    for p in structures:
        new = p.with_name(f'x_{p.name}')
        new.write_bytes(p.read_bytes())
        renamed.append(new)
    second, _, failed = _run(renamed[::-1], tmp_path / 'no_such_usalign', tmp_path / 'c')
    assert not failed
    second = second.rename(index=lambda s: s[2:], columns=lambda s: s[2:])
    pd.testing.assert_frame_equal(first, second.loc[first.index, first.columns])


def test_single_failing_pair_is_isolated(structures, usalign, tmp_path):
    a, b = structures[0].stem, structures[4].stem
    fake = _fake_usalign(tmp_path, usalign,
                         f'case "$*" in */{a}.pdb*/{b}.pdb*|*/{b}.pdb*/{a}.pdb*) exit 1;; esac')
    tm, _, failed = _run(structures, fake, tmp_path / 'c', **PAIRWISE)
    assert [(f[0], f[1]) for f in failed] in ([(a, b)], [(b, a)])
    assert np.isnan(tm.loc[a, b]) and np.isnan(tm.loc[b, a])
    assert tm.isna().sum().sum() == 2


def test_list_mode_skipped_for_mixed_folders(structures, usalign, tmp_path):
    other = tmp_path / 'elsewhere'
    other.mkdir()
    moved = structures[:]
    moved[0] = structures[0].rename(other / structures[0].name)
    fake = _fake_usalign(tmp_path, usalign, 'case "$*" in *-dir*) exit 3;; esac')
    tm, _, failed = _run(moved, fake, tmp_path / 'c', **REAL)
    assert not failed and not tm.isna().any().any()


def test_list_mode_failure_falls_back_to_per_pair(structures, usalign, tmp_path):
    fake = _fake_usalign(tmp_path, usalign, 'case "$*" in *-dir*) exit 1;; esac')
    tm, _, failed = _run(structures, fake, tmp_path / 'c', **REAL)
    ref, _, _ = _run(structures, usalign, tmp_path / 'ref', **PAIRWISE)
    assert not failed
    pd.testing.assert_frame_equal(tm, ref)


def test_end_to_end_with_zip_and_cif(tmp_path, usalign, treecluster, beem):  # pylint: disable=unused-argument
    import zipfile
    from Bio.PDB import MMCIFIO, PDBParser
    from colabalign.pipeline import main
    from conftest import make_structure_set

    src = make_structure_set(tmp_path / 'src')
    indir = tmp_path / 'in'
    indir.mkdir()
    for p in src[:4]:
        (indir / p.name).write_bytes(p.read_bytes())
    with zipfile.ZipFile(indir / 'more.zip', 'w') as z:
        for p in src[4:7]:
            z.write(p, f'batch/{p.name}')
    for p in src[7:]:
        io = MMCIFIO()
        io.set_structure(PDBParser(QUIET=True).get_structure(p.stem, p))
        io.save(str(indir / f'{p.stem}.cif'))

    out = tmp_path / 'out'
    main(['-i', str(indir), '-o', str(out), '-c', '2', '-t', '0.3', '0.6',
          '-u', usalign, '-b', beem, '-m', 'first'])

    results = out / 'results'
    tm = pd.read_csv(results / 'us-align_score_matrix.csv', index_col=0)
    assert len(tm) == len(src)
    assert not tm.isna().any().any()
    assert (results / 'colabalign.tree').stat().st_size > 0
    for threshold in ('0.30', '0.60'):
        assert (results / 'clusters' / threshold).is_dir()
    reps = list(results.rglob('*representatives*'))
    assert reps and all(r.stat().st_size > 0 for r in reps)
    assert list(results.rglob('*_consensus.afasta'))


def test_end_to_end_duplicate_names(tmp_path, usalign, treecluster, beem):  # pylint: disable=unused-argument
    import csv
    import io
    import tarfile
    from colabalign.pipeline import main
    from conftest import make_structure_set

    src = make_structure_set(tmp_path / 'src')
    indir = tmp_path / 'in'
    (indir / 'run1').mkdir(parents=True)
    (indir / 'run2').mkdir()
    # Three different structures that are all called ranked_0.pdb.
    (indir / 'run1' / 'ranked_0.pdb').write_bytes(src[0].read_bytes())
    (indir / 'run2' / 'ranked_0.pdb').write_bytes(src[1].read_bytes())
    with tarfile.open(indir / 'run3.tar.gz', 'w:gz') as t:
        data = src[2].read_bytes()
        info = tarfile.TarInfo('ranked_0.pdb')
        info.size = len(data)
        t.addfile(info, io.BytesIO(data))
    for p in src[3:]:
        (indir / p.name).write_bytes(p.read_bytes())
    originals = {p: p.read_bytes() for p in indir.rglob('*') if p.is_file()}

    out = tmp_path / 'out'
    main(['-i', str(indir / 'run1'), str(indir / 'run2'), str(indir / 'run3.tar.gz'),
          *[str(indir / p.name) for p in src[3:]],
          '-o', str(out), '-c', '2', '-t', '0.3', '-u', usalign, '-b', beem, '-m', 'first'])

    results = out / 'results'
    tm = pd.read_csv(results / 'us-align_score_matrix.csv', index_col=0)
    assert len(tm) == len(src)
    assert {'ranked_0A', 'ranked_0_dup1A', 'ranked_0_dup2A'} <= set(tm.index)
    with open(results / 'duplicates.csv', newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    assert [r['name_in_colabalign'] for r in rows] == [
        'ranked_0A', 'ranked_0_dup1A', 'ranked_0_dup2A']
    assert {r['name_in_colabalign'] for r in rows} <= set(tm.index)
    assert rows[2]['input_file_path'].endswith('run3.tar.gz:ranked_0.pdb')
    assert {p: p.read_bytes() for p in originals} == originals
    assert not (out / '.extracted_inputs').exists()


def test_end_to_end_panel_sizes(tmp_path, usalign, treecluster):  # pylint: disable=unused-argument
    import csv
    from colabalign.pipeline import main
    from conftest import make_structure_set

    src = make_structure_set(tmp_path / 'in', per_family=4)
    out = tmp_path / 'out'
    main(['-i', str(tmp_path / 'in'), '-o', str(out), '-c', '2', '-n', '5', '3-8',
          '-l', '0.5', '-u', usalign])

    results = out / 'results'
    with open(results / 'panel_sizes.csv', newline='', encoding='utf-8') as f:
        sweep = list(csv.DictReader(f))
    assert int(sweep[0]['clusters']) <= len(src) and sweep[-1]['clusters'] == '1'
    requested = [r for r in sweep if r['requested']]
    assert requested and all(r['threshold_used'] for r in requested)
    assert any(3 <= int(r['panel_size']) <= 8 for r in requested)
    # no default 0.20 run when only -n is given, unless the sweep chose it
    with open(results / 'panels.csv', newline='', encoding='utf-8') as f:
        panels = list(csv.DictReader(f))
    assert {p['threshold'] for p in panels} == {r['threshold_used'] for r in requested}
    names = set(pd.read_csv(results / 'us-align_score_matrix.csv', index_col=0).index)
    assert len(names) == len(src)
    for panel in panels:
        with open(results / panel['panel_csv'], newline='', encoding='utf-8') as f:
            rows = list(csv.DictReader(f))
        members = [r['structure'] for r in rows]
        assert len(members) == len(set(members)) == int(panel['panel_size'])
        assert set(members) <= names
        roles = [r['role'] for r in rows]
        assert roles.count('low_coverage') == int(panel['low_coverage_additions'])
        # representatives of multi-member clusters, plus additions: singletons first,
        # then cluster members only once every singleton is in the panel
        tsv = pd.read_csv(results / 'clusters' / f"{panel['threshold']}.tsv", sep='\t')
        singletons = set(tsv.loc[tsv['ClusterNumber'] == -1, 'SequenceName'])
        assert len(singletons) == int(panel['singletons'])
        added = [r['structure'] for r in rows if r['role'] == 'low_coverage']
        assert added[:len(singletons)] == [s for s in added if s in singletons]
        if len(added) > len(singletons):
            assert singletons <= set(added)
        assert int(panel['singletons_in_panel']) == len(singletons & set(added))
        for r in rows:
            if r['structure'] in singletons:
                assert r['role'] == 'low_coverage' and r['cluster'] == '-1'
        assert sum(int(r['structures_represented']) for r in rows) == len(src)
        # every cluster directory has all its members, superposed or copied
        tsv = pd.read_csv(results / 'clusters' / f"{panel['threshold']}.tsv", sep='\t')
        for cluster, group in tsv.groupby('ClusterNumber'):
            pdbs = {p.stem for p in (results / 'clusters' / panel['threshold'] / str(cluster))
                    .glob('*.pdb')}
            assert pdbs == set(group['SequenceName']), (panel['threshold'], cluster)
        fasta = (results / panel['panel_csv']).with_suffix('.fasta').read_text()
        assert fasta.count('>') == len(members)
        assert float(panel['min_best_tm']) <= float(panel['mean_best_tm']) <= 1


def test_end_to_end_default_threshold(tmp_path, usalign, treecluster):  # pylint: disable=unused-argument
    import csv
    from colabalign.pipeline import main
    from conftest import make_structure_set

    make_structure_set(tmp_path / 'in')
    out = tmp_path / 'out'
    main(['-i', str(tmp_path / 'in'), '-o', str(out), '-c', '2', '-u', usalign])
    results = out / 'results'
    assert (results / 'clusters' / '0.20.tsv').is_file()
    assert (results / 'clusters' / '0.20' / 'threshold_0.20_panel.csv').is_file()
    assert (results / 'panels.csv').is_file() and (results / 'panel_sizes.csv').is_file()
    with open(results / 'panels.csv', newline='', encoding='utf-8') as f:
        assert [row['low_coverage_additions'] for row in csv.DictReader(f)] == ['0']  # -l 0 default


def test_no_superpose_gives_same_panels(tmp_path, usalign, treecluster):  # pylint: disable=unused-argument
    from colabalign.pipeline import main
    from conftest import make_structure_set

    make_structure_set(tmp_path / 'in', per_family=4)
    runs = {}
    for name, extra in (('with', []), ('without', ['--no-superpose'])):
        out = tmp_path / name
        main(['-i', str(tmp_path / 'in'), '-o', str(out), '-c', '2', '-n', '4-8',
              '-u', usalign, *extra])
        runs[name] = out / 'results'
    clusters = {name: results / 'clusters' for name, results in runs.items()}
    assert list(clusters['with'].rglob('*.pdb'))
    assert not list(clusters['without'].rglob('*.pdb'))
    for pattern in ('*_cluster_representatives.fasta', '*_panel.csv', '*_panel.fasta', '*.afasta'):
        with_files = sorted(p.relative_to(clusters['with']) for p in clusters['with'].rglob(pattern))
        assert with_files
        assert with_files == sorted(p.relative_to(clusters['without'])
                                    for p in clusters['without'].rglob(pattern))
        for rel in with_files:
            assert (clusters['with'] / rel).read_text() == (clusters['without'] / rel).read_text()
    assert (runs['with'] / 'panels.csv').read_text() == (runs['without'] / 'panels.csv').read_text()
