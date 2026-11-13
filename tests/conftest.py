'''Shared fixtures: synthetic CA-trace structures and external tool discovery.'''

import shutil
from pathlib import Path

import numpy as np
import pytest

RESNAMES = ['ALA', 'GLY', 'LEU', 'SER', 'VAL', 'GLU', 'LYS', 'ASP', 'ILE', 'THR']


def _segment(kind, n, rng):
    '''CA coordinates for an idealised helix, strand or random coil (3.8 A steps).'''
    if kind == 'helix':
        t = np.arange(n) * np.deg2rad(100)
        return np.column_stack([2.3 * np.cos(t), 2.3 * np.sin(t), 1.5 * np.arange(n)])
    if kind == 'strand':
        return np.column_stack([3.3 * np.arange(n), 1.0 * (np.arange(n) % 2), np.zeros(n)])
    steps = rng.normal(size=(n, 3))
    steps = 3.8 * steps / np.linalg.norm(steps, axis=1, keepdims=True)
    return np.cumsum(steps, axis=0)


def make_ca_trace(segments, seed, noise=0.3):
    '''Join segments end to end (each randomly rotated) and add coordinate noise.'''
    rng = np.random.default_rng(seed)
    coords, origin = [], np.zeros(3)
    for kind, n in segments:
        seg = _segment(kind, n, rng)
        q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        seg = (seg - seg[0]) @ q + origin
        coords.append(seg)
        origin = seg[-1] + 3.8 * q[0]
    xyz = np.vstack(coords)
    return xyz + rng.normal(scale=noise, size=xyz.shape)


def write_pdb(path, xyz, chain='A', seed=0, hetatm_calcium=False):
    '''Write a CA-only PDB file (optionally with a trailing HETATM calcium ion).'''
    rng = np.random.default_rng(seed)
    lines = []
    for i, (x, y, z) in enumerate(xyz, start=1):
        res = RESNAMES[rng.integers(len(RESNAMES))]
        lines.append(f'ATOM  {i:5d}  CA  {res} {chain}{i:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C')
    if hetatm_calcium:
        n = len(xyz) + 1
        lines.append(f'HETATM{n:5d} CA    CA {chain}{n:4d}    '
                     f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          CA')
    lines += ['TER', 'END']
    Path(path).write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return Path(path)


# Three loose fold families of different sizes, so clustering has something to find
FAMILIES = {
    'hh': [('helix', 14), ('coil', 4), ('helix', 14)],
    'bab': [('strand', 6), ('coil', 3), ('helix', 12), ('coil', 3), ('strand', 6)],
    'hhh': [('helix', 18), ('coil', 4), ('helix', 18), ('coil', 4), ('helix', 18)],
}


def make_structure_set(directory, per_family=3, extra_residues=(0, 4, 9)):
    '''Write per_family variants of each family (varying length and noise).'''
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for f_idx, (family, segments) in enumerate(FAMILIES.items()):
        for v in range(per_family):
            extra = extra_residues[v % len(extra_residues)]
            segs = segments + ([('coil', extra)] if extra else [])
            seed = 100 * f_idx + v
            xyz = make_ca_trace(segs, seed=seed, noise=0.3 + 0.2 * v)
            paths.append(write_pdb(directory / f'{family}{v}.pdb', xyz, seed=seed))
    return paths


def _tool(*names):
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


@pytest.fixture(scope='session')
def usalign():
    path = _tool('USalign')
    if path is None:
        pytest.skip('USalign not on PATH')
    return path


@pytest.fixture(scope='session')
def treecluster():
    path = _tool('TreeCluster.py')
    if path is None:
        pytest.skip('TreeCluster.py not on PATH')
    return path


@pytest.fixture(scope='session')
def beem():
    path = _tool('BeEM')
    if path is None:
        pytest.skip('BeEM not on PATH')
    return path


@pytest.fixture
def structures(tmp_path):
    '''Nine single-chain synthetic structures in tmp_path/models.'''
    return make_structure_set(tmp_path / 'models')
