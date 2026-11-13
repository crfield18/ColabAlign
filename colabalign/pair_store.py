'''Disk-backed store of per-pair USalign results.

Only TM-scores are kept in memory by the parent process. The per-pair
residue equivalences and superposition matrices live in the on-disk cache
and are loaded on demand, so neither the parent nor any worker ever needs
the results for every pair at once. A PairStore only holds two small
dicts (file hash and canonical rank per structure), so it is cheap to
send to worker processes.
'''

import hashlib
import os
import pickle
from itertools import combinations
from pathlib import Path


def file_hash(path) -> str:
    '''Short sha256 hex digest of a file's contents, used in cache keys.'''
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()[:16]


class PairStore:
    '''
    Per-pair USalign results, cached as one pickle per pair.

    Every pair has a canonical orientation: the structure with the lower
    rank is Structure_1 (the one superposed onto Structure_2). Alignments
    are only ever run, saved and looked up in that orientation, so results
    do not depend on how pairs were batched.
    '''

    def __init__(self, cache_dir: Path, hashes: dict, ranks: dict) -> None:
        '''hashes and ranks map each structure stem to its file hash and canonical rank.'''
        self.cache_dir = Path(cache_dir)
        self.hashes = dict(hashes)
        self.ranks = dict(ranks)

    def ordered(self, stem_a: str, stem_b: str):
        '''Return (first, second, flipped) with first/second in canonical orientation.'''
        if self.ranks[stem_a] <= self.ranks[stem_b]:
            return stem_a, stem_b, False
        return stem_b, stem_a, True

    def cache_file(self, first: str, second: str) -> Path:
        '''Cache path for a pair given in canonical orientation.'''
        return self.cache_dir.joinpath(f'{self.hashes[first]}_{self.hashes[second]}.pkl')

    def load_canonical(self, first: str, second: str):
        '''Load a cached result for a canonically ordered pair, or None if missing/corrupt.'''
        cache_file = self.cache_file(first, second)
        if not cache_file.exists():
            return None
        try:
            with open(cache_file, 'rb') as f:
                return pickle.load(f)
        except Exception:  # pylint: disable=broad-exception-caught
            # any corrupt/unreadable cache entry should just trigger a recompute
            return None

    def load(self, stem_a: str, stem_b: str):
        '''Load a pair in either order. Returns (entry_or_None, flipped).'''
        first, second, flipped = self.ordered(stem_a, stem_b)
        return self.load_canonical(first, second), flipped

    def save(self, first: str, second: str, entry: dict) -> None:
        '''Atomically write a result for a canonically ordered pair.'''
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cache_file = self.cache_file(first, second)
        # Write-then-rename so a reader never sees a half-written file, even
        # if two pairs of identical-content files share a cache key. The
        # .pkl suffix keeps leftovers out of the notebook's results zip.
        tmp_file = cache_file.with_name(f'.{cache_file.stem}.{os.getpid()}.tmp.pkl')
        with open(tmp_file, 'wb') as f:
            pickle.dump(entry, f)
        os.replace(tmp_file, cache_file)

    def transform(self, source: str, target: str):
        '''3x4 matrix superposing source onto target, or None if the pair is missing.'''
        entry, flipped = self.load(source, target)
        if entry is None:
            return None
        return entry['transform_mx_reverse' if flipped else 'transform_mx_forward']

    def sub_library(self, members) -> dict:
        '''
        Residue-equivalence library restricted to pairs within members, keyed
        by canonically ordered (first, second) stems. This is all a cluster's
        MSA needs, so cluster jobs never load anything outside their cluster.
        '''
        library = {}
        for stem_a, stem_b in combinations(members, 2):
            first, second, _ = self.ordered(stem_a, stem_b)
            entry = self.load_canonical(first, second)
            if entry is None:
                continue
            library[(first, second)] = {
                'res_a': entry['res_a'],
                'res_b': entry['res_b'],
                'weight': entry['weight'],
            }
        return library
