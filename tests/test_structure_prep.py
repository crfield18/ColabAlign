'''Tests for sequence extraction and input discovery (including .zip inputs).'''

import zipfile

from colabalign.common import count_residues, extract_sequence
from colabalign.structure_prep import discover_input_models
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
    assert {'ranked_0', 'single', 'model', 'run1_ranked_0', 'run2_ranked_0'} == set(stems)
    assert all(m.exists() for m in models)
    assert not any('__MACOSX' in str(m) or m.name.startswith('.') for m in models)


def test_zip_given_directly(tmp_path):
    z = _zip(tmp_path / 'structures.zip', {'a.pdb': 'END\n', 'b/c.cif': 'data_c\n'})
    models = discover_input_models([z], tmp_path / 'extract')
    assert sorted(m.name for m in models) == ['a.pdb', 'c.cif']
