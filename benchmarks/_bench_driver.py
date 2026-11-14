'''Helper run by benchmark_versions.py: prepare structures or run only the
alignment step of a given ColabAlign source tree.

    python _bench_driver.py prepare --src SRC --inputs DIR --models DIR --beem BEEM --mode first --cores N
    python _bench_driver.py align   --src SRC --models DIR --out DIR --usalign USALIGN --cores N
    python _bench_driver.py rawtable --table FILE --models DIR --out DIR

It works with both the flat layout (<= 1.1.x, modules at the repo root) and
the colabalign/ package (>= 1.2.0).
'''

import argparse
import shutil
import sys
import tempfile
import time
from pathlib import Path


def _import(src: Path, module: str):
    '''Import a ColabAlign module from src, whichever layout it uses.'''
    sys.path.insert(0, str(src))
    if src.joinpath('colabalign', '__init__.py').exists():
        return __import__(f'colabalign.{module}', fromlist=[module])
    return __import__(module)


def prepare(args):
    '''Run structure preparation once, leaving single-chain PDBs in --models.'''
    prep = _import(args.src, 'structure_prep')
    models = Path(args.models)
    tmp = models.with_name(models.name + '.tmp')
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    with tempfile.TemporaryDirectory() as extract:
        try:
            inputs = prep.discover_input_models([Path(args.inputs)], Path(extract))
        except TypeError:  # older signature without an extraction folder
            inputs = prep.discover_input_models([Path(args.inputs)])
        prep.prepare_models(inputs, tmp, Path(args.beem), args.mode, args.cores)
    tmp.rename(models)


def align(args):
    '''Time only the all-vs-all alignment step and write the TM-score matrix.'''
    runner = _import(args.src, 'usalign_runner')
    models = Path(args.models)
    model_list = sorted(models.glob('*.pdb'))
    out = Path(args.out)
    results = out / 'results'
    results.mkdir(parents=True, exist_ok=True)
    cache = results / 'usalign_cache'

    start = time.perf_counter()
    try:  # >= 1.2.0
        result = runner.run_pairwise_alignments(model_list, args.usalign, args.cores, cache)
    except TypeError:  # <= 1.1.x also takes models_path
        result = runner.run_pairwise_alignments(model_list, models, args.usalign, args.cores, cache)
    elapsed = time.perf_counter() - start
    print(f'Pairwise alignment elapsed time:\t{elapsed:.2f} s')
    result[0].to_csv(results / 'us-align_score_matrix.csv', index=True)


def rawtable(args):
    '''Convert USalign -outfmt 2 output into a TM-score matrix CSV (not timed).'''
    import numpy as np  # pylint: disable=import-outside-toplevel
    import pandas as pd  # pylint: disable=import-outside-toplevel

    names = sorted(p.stem for p in Path(args.models).glob('*.pdb'))
    index = {n: i for i, n in enumerate(names)}
    matrix = np.full((len(names), len(names)), np.nan)
    np.fill_diagonal(matrix, 1.0)

    def stem(field):
        return Path(field.split(':', 1)[0]).stem

    with open(args.table, encoding='utf-8', errors='ignore') as f:
        for line in f:
            parts = line.split('\t')
            if len(parts) < 4 or line.startswith('#'):
                continue
            try:
                i, j = index[stem(parts[0])], index[stem(parts[1])]
                tm = max(float(parts[2]), float(parts[3]))
            except (KeyError, ValueError):
                continue
            matrix[i, j] = matrix[j, i] = tm
    results = Path(args.out) / 'results'
    results.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(matrix, index=names, columns=names).to_csv(results / 'us-align_score_matrix.csv')


def main():
    '''Dispatch the sub-command.'''
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    s = sub.add_parser('prepare')
    s.add_argument('--src', type=Path, required=True)
    s.add_argument('--inputs', required=True)
    s.add_argument('--models', required=True)
    s.add_argument('--beem', required=True)
    s.add_argument('--mode', default='first')
    s.add_argument('--cores', type=int, default=1)
    s = sub.add_parser('align')
    s.add_argument('--src', type=Path, required=True)
    s.add_argument('--models', required=True)
    s.add_argument('--out', required=True)
    s.add_argument('--usalign', required=True)
    s.add_argument('--cores', type=int, default=1)
    s = sub.add_parser('rawtable')
    s.add_argument('--table', required=True)
    s.add_argument('--models', required=True)
    s.add_argument('--out', required=True)
    args = p.parse_args()
    {'prepare': prepare, 'align': align, 'rawtable': rawtable}[args.cmd](args)


if __name__ == '__main__':
    main()
