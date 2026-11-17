'''Tests for sequence extraction and input discovery (including .zip inputs).'''

import zipfile

import pytest

from colabalign.common import count_residues, extract_sequence
import csv

from colabalign.structure_prep import (
    collect_input_models, discover_input_models, write_duplicates_csv, DUPLICATES_CSV_COLUMNS,
)
from conftest import make_ca_trace, write_pdb


def test_hetatm_calcium_is_not_a_residue(tmp_path):
    xyz = make_ca_trace([('helix', 12)], seed=0)
    plain = write_pdb(tmp_path / 'plain.pdb', xyz, seed=1)
    with_ca = write_pdb(tmp_path / 'ca.pdb', xyz, seed=1, hetatm_calcium=True)
    assert extract_sequence(with_ca) == extract_sequence(plain)
    assert len(extract_sequence(plain)) == 12 and 'X' not in extract_sequence(plain)
    assert count_residues(with_ca) == 12


def test_mse_hetatm_is_kept(tmp_path):
    path = tmp_path / 'mse.pdb'
    path.write_text(
        'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n'
        'HETATM    2  CA  MSE A   2       3.800   0.000   0.000  1.00  0.00           C\n'
        'END\n')
    assert len(extract_sequence(path)) == 2


def test_count_residues_minimum_one(tmp_path):
    path = tmp_path / 'empty.pdb'
    path.write_text('END\n')
    assert count_residues(path) == 1


def _zip(path, members):
    with zipfile.ZipFile(path, 'w') as z:
        for name, data in members.items():
            z.writestr(name, data)
    return path


def test_discover_mixed_inputs_and_zips(tmp_path):
    pdb = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n'
    indir = tmp_path / 'in'
    indir.mkdir()
    (indir / 'ranked_0.pdb').write_text(pdb)
    (indir / 'notes.txt').write_text('x')
    (indir / '._ranked_0.pdb').write_text('resource fork')
    _zip(indir / 'runs.ZIP', {
        'run1/ranked_0.pdb': pdb,           # clashes with the loose file
        'run2/ranked_0.pdb': pdb,           # clashes with run1 as well
        'run2/model.CIF': 'data_x\n',
        '__MACOSX/run1/._ranked_0.pdb': 'junk',
        '.hidden/x.pdb': pdb,
        'readme.md': 'skip',
    })
    (indir / 'broken.zip').write_bytes(b'not a zip')
    loose = tmp_path / 'single.pdb'
    loose.write_text(pdb)

    models = discover_input_models([indir, loose, tmp_path / 'missing'], tmp_path / 'extract')
    stems = sorted(m.stem.lower() for m in models)
    assert len(stems) == len(set(stems)) == 5
    assert {'ranked_0', 'single', 'model', 'ranked_0_dup1', 'ranked_0_dup2'} == set(stems)
    assert all(m.exists() for m in models)
    assert not any('__MACOSX' in str(m) or m.name.startswith('.') for m in models)


def test_zip_given_directly(tmp_path):
    z = _zip(tmp_path / 'structures.zip', {'a.pdb': 'END\n', 'b/c.cif': 'data_c\n'})
    models = discover_input_models([z], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['a.pdb', 'c.cif']


def _tar(path, members, mode='w:gz'):
    import io
    import tarfile
    with tarfile.open(path, mode) as t:
        for name, data in members.items():
            raw = data.encode()
            info = tarfile.TarInfo(name)
            info.size = len(raw)
            t.addfile(info, io.BytesIO(raw))
    return path


@pytest.mark.parametrize('name,mode', [
    ('structures.tar.gz', 'w:gz'), ('structures.TGZ', 'w:gz'),
    ('structures.tar.bz2', 'w:bz2'), ('structures.tar.xz', 'w:xz'),
    ('structures.tar', 'w'),
])
def test_tar_given_directly(tmp_path, name, mode):
    t = _tar(tmp_path / name, {'a.pdb': 'END\n', 'b/c.cif': 'data_c\n',
                               '._a.pdb': 'junk', 'notes.txt': 'skip'}, mode)
    models = discover_input_models([t], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['a.pdb', 'c.cif']
    # extracted under the full archive stem, not 'structures.tar'
    assert all((tmp_path / 'extract' / 'structures') in m.parents for m in models)


def test_tar_in_directory_dup_suffixes_clashes_and_skips_broken(tmp_path):
    pdb = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n'
    indir = tmp_path / 'in'
    indir.mkdir()
    (indir / 'ranked_0.pdb').write_text(pdb)
    _tar(indir / 'runs.tar.gz', {'./run1/ranked_0.pdb': pdb, 'run2/ranked_0.pdb': pdb,
                                 '__MACOSX/run1/._ranked_0.pdb': 'junk'})
    (indir / 'broken.tar.gz').write_bytes(b'not a tarball')
    models = discover_input_models([indir], tmp_path / 'extract')
    assert sorted(m.stem for m in models) == ['ranked_0', 'ranked_0_dup1', 'ranked_0_dup2']
    assert all(m.exists() for m in models)


def test_tar_rejects_path_traversal_and_links(tmp_path):
    import io
    import tarfile
    path = tmp_path / 'evil.tar.gz'
    with tarfile.open(path, 'w:gz') as t:
        for name in ('../escape.pdb', '/abs.pdb', 'ok.pdb'):
            info = tarfile.TarInfo(name)
            info.size = 4
            t.addfile(info, io.BytesIO(b'END\n'))
        link = tarfile.TarInfo('link.pdb')
        link.type, link.linkname = tarfile.SYMTYPE, '/etc/passwd'
        t.addfile(link)
    models = discover_input_models([path], tmp_path / 'extract')
    assert [m.name for m in models] == ['ok.pdb']
    assert not (tmp_path / 'escape.pdb').exists()
    assert not any(p.is_symlink() for p in (tmp_path / 'extract').rglob('*'))


PDB_A = 'ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n'
PDB_B = 'ATOM      1  CA  ALA A   1       1.234   0.000   0.000  1.00  0.00           C\nEND\n'


def test_loose_duplicates_get_dup_copies_and_originals_untouched(tmp_path):
    # Same byte size, different coordinates: the old size check treated these as one file.
    assert len(PDB_A) == len(PDB_B)
    for sub, text in (('A', PDB_A), ('B', PDB_B), ('C', PDB_A)):
        (tmp_path / sub).mkdir()
        (tmp_path / sub / 'ranked_0.pdb').write_text(text)
    before = {p: p.read_text() for p in tmp_path.rglob('*.pdb')}

    models, rows = collect_input_models(
        [tmp_path / 'A', tmp_path / 'B', tmp_path / 'C'], tmp_path / 'extract')

    assert sorted(m.name for m in models) == ['ranked_0.pdb', 'ranked_0_dup1.pdb',
                                              'ranked_0_dup2.pdb']
    # originals unchanged, copies hold the right contents
    assert {p: p.read_text() for p in before} == before
    by_name = {m.name: m for m in models}
    assert by_name['ranked_0.pdb'] == tmp_path / 'A' / 'ranked_0.pdb'
    assert by_name['ranked_0_dup1.pdb'].read_text() == PDB_B
    assert by_name['ranked_0_dup2.pdb'].read_text() == PDB_A
    assert all(tmp_path / 'extract' in by_name[n].parents
               for n in ('ranked_0_dup1.pdb', 'ranked_0_dup2.pdb'))

    assert [(r['input_file_path'], r['name_in_colabalign'], r['identical_to']) for r in rows] == [
        (str((tmp_path / 'A' / 'ranked_0.pdb').resolve()), 'ranked_0.pdb', ''),
        (str((tmp_path / 'B' / 'ranked_0.pdb').resolve()), 'ranked_0_dup1.pdb', ''),
        (str((tmp_path / 'C' / 'ranked_0.pdb').resolve()), 'ranked_0_dup2.pdb', 'ranked_0.pdb'),
    ]
    assert all(r['input_file_name'] == 'ranked_0.pdb' for r in rows)
    assert rows[0]['sha256'] == rows[2]['sha256'] != rows[1]['sha256']


def test_archive_duplicates_report_member_paths(tmp_path):
    (tmp_path / 'ranked_0.pdb').write_text(PDB_A)
    archive = _tar(tmp_path / 'runs.tar.gz', {'run1/ranked_0.pdb': PDB_B, 'run1/other.pdb': PDB_A})
    models, rows = collect_input_models([tmp_path / 'ranked_0.pdb', archive], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['other.pdb', 'ranked_0.pdb', 'ranked_0_dup1.pdb']
    assert [r['input_file_path'] for r in rows] == [
        str((tmp_path / 'ranked_0.pdb').resolve()), f'{archive.resolve()}:run1/ranked_0.pdb']
    assert 'other.pdb' not in {r['input_file_name'] for r in rows}  # no clash, not listed


def test_stem_clash_across_extensions_and_case(tmp_path):
    (tmp_path / 'model.pdb').write_text(PDB_A)
    (tmp_path / 'MODEL.cif').write_text('data_x\n')
    models, rows = collect_input_models([tmp_path], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['MODEL.cif', 'model_dup1.pdb']
    assert len(rows) == 2


def test_dup_suffix_skips_names_already_taken(tmp_path):
    for sub in 'AB':
        (tmp_path / sub).mkdir()
        (tmp_path / sub / 'x.pdb').write_text(PDB_A)
    (tmp_path / 'A' / 'x_dup1.pdb').write_text(PDB_B)  # a real input already called x_dup1
    models, _ = collect_input_models([tmp_path / 'A', tmp_path / 'B'], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['x.pdb', 'x_dup1.pdb', 'x_dup2.pdb']


def test_same_file_given_twice_is_not_a_duplicate(tmp_path):
    (tmp_path / 'a.pdb').write_text(PDB_A)
    models, rows = collect_input_models([tmp_path, tmp_path / 'a.pdb'], tmp_path / 'extract')
    assert [m.name for m in models] == ['a.pdb'] and rows == []


def test_write_duplicates_csv(tmp_path):
    for sub in 'AB':
        (tmp_path / sub).mkdir()
        (tmp_path / sub / 'x.pdb').write_text(PDB_A)
    _, rows = collect_input_models([tmp_path / 'A', tmp_path / 'B'], tmp_path / 'extract')
    write_duplicates_csv(rows, tmp_path / 'duplicates.csv')
    with open(tmp_path / 'duplicates.csv', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        assert tuple(reader.fieldnames) == DUPLICATES_CSV_COLUMNS
        assert [r['name_in_colabalign'] for r in reader] == ['x.pdb', 'x_dup1.pdb']

    write_duplicates_csv([], tmp_path / 'empty.csv')
    assert (tmp_path / 'empty.csv').read_text().strip() == ','.join(DUPLICATES_CSV_COLUMNS)


@pytest.mark.parametrize('mode', ['first', 'all'])
def test_same_size_duplicates_survive_preparation(tmp_path, mode):
    from colabalign.structure_prep import prepare_models
    for sub, text in (('A', PDB_A), ('B', PDB_B)):
        (tmp_path / sub).mkdir()
        (tmp_path / sub / 'ranked_0.pdb').write_text(text)
    models = discover_input_models([tmp_path / 'A', tmp_path / 'B'], tmp_path / 'extract')
    out = tmp_path / 'models'
    out.mkdir()
    final = prepare_models(models, out, tmp_path / 'no_beem', mode, 1)
    assert sorted(p.name for p in final) == ['ranked_0A.pdb', 'ranked_0_dup1A.pdb']


def test_safe_copy_destination_uses_content_not_size(tmp_path):
    from colabalign.structure_prep import _safe_copy_destination
    dst, same, diff = tmp_path / 'm.pdb', tmp_path / 's.pdb', tmp_path / 'd.pdb'
    dst.write_text(PDB_A)
    same.write_text(PDB_A)
    diff.write_text(PDB_B)
    assert _safe_copy_destination(dst, same) == dst
    assert _safe_copy_destination(dst, diff) == tmp_path / 'm_dup1.pdb'
