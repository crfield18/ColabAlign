'''Integration tests against the real external tools (skipped when not installed).'''

import stat

import numpy as np
import pandas as pd
import pytest

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
    for key, entry in la.items():
        for field, value in entry.items():
            assert np.array_equal(value, lb[key][field])


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
