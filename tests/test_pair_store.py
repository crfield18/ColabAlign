'''Tests for the on-disk PairStore and transform handling.'''
# pylint: disable=redefined-outer-name

import numpy as np
import pytest

from colabalign.pair_store import PairStore, file_hash
from colabalign.usalign_runner import reverse_transformation_matrix


def _random_transform(seed):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return np.column_stack([rng.normal(size=3) * 10, q])  # [t | U], USalign layout


def _apply(mx, xyz):
    return xyz @ mx[:, 1:].T + mx[:, 0]


def test_reverse_transformation_round_trip():
    mx = _random_transform(1)
    xyz = np.random.default_rng(2).normal(size=(20, 3)) * 10
    back = _apply(reverse_transformation_matrix(mx), _apply(mx, xyz))
    assert np.allclose(back, xyz)


def test_file_hash_depends_on_content_only(tmp_path):
    a, b, c = tmp_path / 'a.pdb', tmp_path / 'b.pdb', tmp_path / 'c.pdb'
    a.write_text('ATOM 1\n')
    b.write_text('ATOM 1\n')
    c.write_text('ATOM 2\n')
    assert file_hash(a) == file_hash(b) != file_hash(c)
    assert len(file_hash(a)) == 16


@pytest.fixture
def store(tmp_path):
    stems = ['long', 'mid', 'short']
    return PairStore(tmp_path / 'cache',
                     hashes={s: f'h{s}' for s in stems},
                     ranks={s: i for i, s in enumerate(stems)})


def _entry(seed):
    fwd = _random_transform(seed)
    return {'tm1': 0.5, 'tm2': 0.6, 'res_a': [0, 1, 2], 'res_b': [3, 4, 5],
            'weight': [1.0, 0.5, 0.2], 'transform_mx_forward': fwd,
            'transform_mx_reverse': reverse_transformation_matrix(fwd)}


def test_orientation(store):
    assert store.ordered('long', 'short') == ('long', 'short', False)
    assert store.ordered('short', 'long') == ('long', 'short', True)


def test_save_load_and_transform_direction(store):
    entry = _entry(3)
    store.save('long', 'short', entry)
    assert store.cache_file('long', 'short').exists()
    assert not list(store.cache_dir.glob('.*tmp*'))          # atomic write left nothing behind
    loaded, flipped = store.load('short', 'long')
    assert flipped and loaded['res_a'] == [0, 1, 2]
    assert np.allclose(store.transform('long', 'short'), entry['transform_mx_forward'])
    assert np.allclose(store.transform('short', 'long'), entry['transform_mx_reverse'])
    assert store.transform('long', 'mid') is None


def test_corrupt_cache_entry_is_treated_as_missing(store):
    store.cache_dir.mkdir(parents=True)
    store.cache_file('long', 'mid').write_bytes(b'not a pickle')
    assert store.load_canonical('long', 'mid') is None


def test_sub_library_only_loads_members(store):
    store.save('long', 'mid', _entry(1))
    store.save('long', 'short', _entry(2))
    store.save('mid', 'short', _entry(3))
    lib = store.sub_library(['short', 'long'])   # any member order
    assert set(lib) == {('long', 'short')}
    assert set(lib[('long', 'short')]) == {'res_a', 'res_b', 'weight'}
    assert set(store.sub_library(['long', 'mid', 'short'])) == {
        ('long', 'mid'), ('long', 'short'), ('mid', 'short')}
