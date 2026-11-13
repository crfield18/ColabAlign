#!/usr/bin/env python3
'''Benchmark ColabAlign versions on nested random subsets of a structure set.

Two stages can be benchmarked (--stages, default both):

  full   the whole pipeline (structure preparation, alignment, tree,
         clustering, representatives) for the old and new versions.
  align  only the all-vs-all alignment step. Structures for each subset
         are prepared once (with the new version's code), then the old
         version, the new version and a single plain `USalign -dir`
         process all align the same chain files.

Subsets are nested (the 50-structure set contains the 20-structure set),
and the run order of the contenders rotates at each size so that drift
(thermal throttling, other load) does not favour any of them.

Per run it records wall time, CPU time, the phase timings ColabAlign
prints, peak memory of the main process and of the whole process tree,
the cache size on disk, and the number and total length of chains. For
each stage and size it also compares every contender's TM-score matrix
with the new version's.

Results are appended to <workdir>/results.csv as runs finish, so an
interrupted benchmark resumes when the same command is run again.

Example (macOS: caffeinate stops the machine sleeping mid-run):

    caffeinate -i python3 benchmarks/benchmark_versions.py \\
        --structures ~/data/pdbs --workdir ~/colabalign-bench

The python3 used must have ColabAlign's dependencies installed (for example
the colabalign conda environment), and USalign, BeEM and TreeCluster.py
must be on PATH or given with --usalign / --beem.
'''

import argparse
import csv
import io
import os
import random
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import time
from pathlib import Path

DEFAULT_SIZES = [10, 20, 50, 100, 200, 300, 500, 750, 1000, 1200]
STRUCTURE_SUFFIXES = {'.pdb', '.cif'}
DRIVER = Path(__file__).resolve().with_name('_bench_driver.py')
PHASES = {
    'Structure preparation': 'prep_s',
    'Pairwise alignment': 'align_s',
    'Dendrogram calculation': 'tree_s',
    'Dendrogram clustering': 'cluster_s',
    'Representative calculation': 'reps_s',
    'Total': 'total_reported_s',
}
FIELDS = ['stage', 'version', 'ref', 'size', 'repeat', 'status', 'chains', 'residues',
          'wall_s', 'cpu_s', *PHASES.values(), 'peak_main_rss_mb', 'peak_tree_rss_mb',
          'peak_processes', 'cache_mb', 'cache_files', 'log']
STAGE_CONTENDERS = {'full': ['old', 'new'], 'align': ['old', 'new', 'usalign_dir']}
LABELS = {'old': 'baseline', 'new': 'new', 'usalign_dir': 'USalign -dir (1 process)'}
COLOURS = {'old': 'tab:blue', 'new': 'tab:orange', 'usalign_dir': 'tab:green'}


def parse_args():
    '''Command-line options.'''
    p = argparse.ArgumentParser(description=__doc__.split('\n\n', maxsplit=1)[0])
    p.add_argument('--structures', type=Path, required=True,
                   help='Directory of .pdb/.cif files to sample from.')
    p.add_argument('--workdir', type=Path, default=Path('colabalign-benchmark'),
                   help='Where sources, subsets, logs and results go (default ./colabalign-benchmark).')
    p.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[1],
                   help='ColabAlign git checkout (default: the one containing this script).')
    p.add_argument('--old-ref', default='3077acf',
                   help='Baseline git ref (default 3077acf, i.e. v1.1.1 plus the two one-line '
                        'fixes it needs to run; the v1.1.1 tag itself has a syntax error).')
    p.add_argument('--new-ref', default='HEAD',
                   help="Git ref to compare (default HEAD). Use 'WORKTREE' for uncommitted changes.")
    p.add_argument('--stages', nargs='+', choices=['full', 'align'], default=['full', 'align'],
                   help='What to benchmark (default both).')
    p.add_argument('--contenders', nargs='+', choices=['old', 'new', 'usalign_dir'],
                   default=['old', 'new', 'usalign_dir'],
                   help='Restrict which contenders run (usalign_dir only runs in the align stage).')
    p.add_argument('--raw-opts', default='-het 2 -outfmt 2',
                   help="Options for the plain USalign -dir run (default '-het 2 -outfmt 2': same "
                        "residues as ColabAlign, scores only). Use '-het 2 -do -m -' to make "
                        'USalign produce everything ColabAlign asks it for.')
    p.add_argument('--sizes', type=int, nargs='+', default=DEFAULT_SIZES,
                   help=f'Subset sizes (default {" ".join(map(str, DEFAULT_SIZES))}). '
                        'Sizes larger than the set are replaced by the full set.')
    p.add_argument('--repeats', type=int, default=1, help='Runs per contender per size (default 1).')
    p.add_argument('--cores', type=int, default=os.cpu_count(), help='Cores per run (default all).')
    p.add_argument('--threshold', type=float, nargs='+', default=[0.20])
    p.add_argument('--mode', choices=['first', 'all'], default='first')
    p.add_argument('--seed', type=int, default=1, help='Subset sampling seed (default 1).')
    p.add_argument('--timeout', type=float, default=6.0,
                   help='Per-run limit in hours (default 6). A contender that times out or fails '
                        'is not run at larger sizes in that stage.')
    p.add_argument('--max-rss-gb', type=float, default=None,
                   help='Kill a run whose process tree exceeds this much memory.')
    p.add_argument('--usalign', default='USalign')
    p.add_argument('--beem', default='BeEM')
    p.add_argument('--python', default=sys.executable,
                   help='Interpreter used to run ColabAlign (default: this one).')
    p.add_argument('--start-method', choices=['auto', 'default', 'spawn', 'fork'], default='auto',
                   help="multiprocessing start method for ColabAlign runs. 'auto' (default) only "
                        "switches to spawn when the interpreter defaults to forkserver (Linux, "
                        'Python >= 3.14), whose workers are not counted in CPU time.')
    p.add_argument('--keep-outputs', action='store_true',
                   help="Keep each run's output folder (the per-pair cache can be many GB).")
    return p.parse_args()


# Sources

def export_source(repo: Path, ref: str, dest: Path) -> Path:
    '''Export a git ref (or the working tree) into dest and return it.'''
    if dest.exists():
        return dest
    if ref == 'WORKTREE':
        shutil.copytree(repo, dest, ignore=shutil.ignore_patterns('.git', '__pycache__', '*.egg-info'))
        return dest
    sha = subprocess.run(['git', '-C', str(repo), 'rev-parse', '--verify', f'{ref}^{{commit}}'],
                         capture_output=True, text=True, check=True).stdout.strip()
    archive = subprocess.run(['git', '-C', str(repo), 'archive', '--format=tar', sha],
                             capture_output=True, check=True).stdout
    tmp = dest.with_name(dest.name + '.tmp')
    shutil.rmtree(tmp, ignore_errors=True)
    tmp.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(archive)) as tar:
        tar.extractall(tmp)  # trusted: our own repository
    tmp.rename(dest)
    return dest


def pipeline_command(src: Path, python: str):
    '''Command prefix and extra environment to run the full pipeline from a source tree.'''
    if src.joinpath('colabalign', '__main__.py').exists():  # packaged (>= 1.2.0)
        return [python, '-u', '-m', 'colabalign'], {'PYTHONPATH': str(src)}
    return [python, '-u', str(src / 'colabalign.py')], {}  # flat layout (<= 1.1.x)


# Subsets

def make_subsets(structures: Path, sizes, seed: int, workdir: Path):
    '''Nested random subsets as folders of symlinks. Returns {size: folder}.'''
    files = sorted(f for f in structures.iterdir()
                   if f.is_file() and f.suffix.lower() in STRUCTURE_SUFFIXES
                   and not f.name.startswith('.'))
    if not files:
        sys.exit(f'No .pdb/.cif files found in {structures}')
    random.Random(seed).shuffle(files)
    sizes = sorted({min(s, len(files)) for s in sizes if s >= 2})
    subsets = {}
    for size in sizes:
        folder = workdir / 'inputs' / f'n{size:05d}'
        if not folder.exists():
            tmp = folder.with_name(folder.name + '.tmp')
            shutil.rmtree(tmp, ignore_errors=True)
            tmp.mkdir(parents=True)
            for f in files[:size]:
                tmp.joinpath(f.name).symlink_to(f.resolve())
            tmp.rename(folder)
        subsets[size] = folder
    print(f'{len(files)} structures available; sizes: {" ".join(map(str, sizes))}')
    return subsets


# Measurement

def process_table():
    '''{pid: (ppid, rss_kb)} for every process (ps works on macOS and Linux).'''
    out = subprocess.run(['ps', '-A', '-o', 'pid=,ppid=,rss='],
                         capture_output=True, text=True, check=False).stdout
    table = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) == 3 and all(x.isdigit() for x in parts):
            table[int(parts[0])] = (int(parts[1]), int(parts[2]))
    return table


def tree_usage(root: int):
    '''(main RSS MB, summed tree RSS MB, process count) for root and its descendants.
    Summed RSS counts pages shared between forked workers more than once.'''
    table = process_table()
    children = {}
    for pid, (ppid, _) in table.items():
        children.setdefault(ppid, []).append(pid)
    stack, total, count = [root], 0, 0
    while stack:
        pid = stack.pop()
        if pid in table:
            total += table[pid][1]
            count += 1
        stack.extend(children.get(pid, []))
    main_kb = table.get(root, (0, 0))[1]
    return main_kb / 1024, total / 1024, count


def children_cpu():
    '''User + system CPU seconds of all reaped descendants so far.'''
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def dir_usage(path: Path):
    '''(MB, file count) of a folder tree.'''
    total = count = 0
    for root, _, names in os.walk(path):
        for name in names:
            try:
                total += os.lstat(os.path.join(root, name)).st_size
                count += 1
            except OSError:
                pass
    return total / 1e6, count


def chain_stats(models: Path):
    '''(chains, total residues) of the prepared single-chain models.'''
    chains = residues = 0
    for f in models.glob('*.pdb'):
        chains += 1
        seen = set()
        with open(f, encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                    seen.add((line[21], line[22:27]))
        residues += len(seen)
    return chains, residues


def parse_phases(log_text: str):
    '''Phase timings from ColabAlign's "<phase> elapsed time:\\t<s> s" lines.'''
    phases = {}
    for line in log_text.replace('\r', '\n').splitlines():
        if 'elapsed time:' not in line:
            continue
        label, _, value = line.partition('elapsed time:')
        key = PHASES.get(label.strip())
        try:
            if key:
                phases[key] = round(float(value.strip().split()[0]), 2)
        except (ValueError, IndexError):
            pass
    return phases




# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def run_measured(cmd, env_extra, log_path, args, stdout_path=None):
    '''Run a command, polling memory of its process tree. Returns measurements.'''
    env = {**os.environ, **env_extra, 'PYTHONUNBUFFERED': '1', 'TQDM_MININTERVAL': '10'}
    limit_s = args.timeout * 3600
    peak_main = peak_tree = peak_procs = 0.0
    status = 'ok'
    cpu0, t0 = children_cpu(), time.perf_counter()
    with open(log_path, 'w', encoding='utf-8') as log:
        log.write(shlex.join(map(str, cmd)) + '\n\n')
        log.flush()
        stdout = open(stdout_path, 'w', encoding='utf-8') if stdout_path else log  # pylint: disable=consider-using-with
        try:
            proc = subprocess.Popen(cmd, stdout=stdout, stderr=log, env=env,  # pylint: disable=consider-using-with
                                    start_new_session=True)
            try:
                while proc.poll() is None:
                    main_mb, tree, procs = tree_usage(proc.pid)
                    peak_main, peak_tree = max(peak_main, main_mb), max(peak_tree, tree)
                    peak_procs = max(peak_procs, procs)
                    if time.perf_counter() - t0 > limit_s:
                        status = 'timeout'
                    elif args.max_rss_gb and tree > args.max_rss_gb * 1024:
                        status = 'memory_limit'
                    if status != 'ok':
                        os.killpg(proc.pid, signal.SIGKILL)
                        break
                    time.sleep(0.5)
            except KeyboardInterrupt:
                os.killpg(proc.pid, signal.SIGKILL)
                proc.wait()
                raise
            proc.wait()
        finally:
            if stdout_path:
                stdout.close()
    wall = time.perf_counter() - t0
    if status == 'ok' and proc.returncode != 0:
        status = f'exit_{proc.returncode}'
    row = {'status': status, 'wall_s': round(wall, 2), 'cpu_s': round(children_cpu() - cpu0, 2),
           'peak_main_rss_mb': round(peak_main, 1), 'peak_tree_rss_mb': round(peak_tree, 1),
           'peak_processes': int(peak_procs), 'log': str(log_path)}
    row.update(parse_phases(log_path.read_text(encoding='utf-8', errors='ignore')))
    return row


class Bench:
    '''Holds sources, inputs and settings, and runs one contender at a time.'''

    def __init__(self, args, workdir: Path):
        self.args = args
        self.workdir = workdir
        self.sources = {}
        for version, ref in (('old', args.old_ref), ('new', args.new_ref)):
            src = export_source(args.repo, ref, workdir / 'src' / f'{version}_{ref.replace("/", "_")}')
            self.sources[version] = (ref, src)
            print(f'{version}: {ref} -> {src}')
        self.refs = {'old': args.old_ref, 'new': args.new_ref, 'usalign_dir': 'USalign'}
        self.pythonpath = self._start_method_hook()

    def _start_method_hook(self):
        '''Folder with a sitecustomize.py that sets the multiprocessing start method, or None.'''
        method = self.args.start_method
        if method == 'auto':
            default = subprocess.run(
                [self.args.python, '-c', 'import multiprocessing as m; print(m.get_start_method())'],
                capture_output=True, text=True, check=True).stdout.strip()
            method = 'spawn' if default == 'forkserver' else 'default'
            if method == 'spawn':
                print('Interpreter defaults to forkserver; using spawn so worker CPU time is counted.')
        if method == 'default':
            return None
        hook = self.workdir / 'start_method_hook'
        hook.mkdir(exist_ok=True)
        hook.joinpath('sitecustomize.py').write_text(
            'import multiprocessing\n'
            f'multiprocessing.set_start_method({method!r}, force=True)\n', encoding='utf-8')
        return hook

    def python_env(self, env):
        '''Add the start-method hook to a run's PYTHONPATH.'''
        if self.pythonpath is None:
            return env
        paths = [str(self.pythonpath)] + ([env['PYTHONPATH']] if 'PYTHONPATH' in env else [])
        return {**env, 'PYTHONPATH': os.pathsep.join(paths)}

    def prepared_models(self, size, inputs):
        '''Single-chain models for a subset, prepared once with the new version.'''
        models = self.workdir / 'prepared' / f'n{size:05d}'
        if not models.exists():
            log = self.workdir / 'logs' / f'prepare_n{size:05d}.log'
            cmd = [self.args.python, str(DRIVER), 'prepare', '--src', str(self.sources['new'][1]),
                   '--inputs', str(inputs), '--models', str(models), '--beem', self.args.beem,
                   '--mode', self.args.mode, '--cores', str(self.args.cores)]
            with open(log, 'w', encoding='utf-8') as f:
                subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT, check=True)
        return models

    def run(self, stage, version, size, inputs, out, log):
        '''Run one contender for one stage and subset. Returns a result row.'''
        args = self.args
        shutil.rmtree(out, ignore_errors=True)
        out.mkdir(parents=True)
        if stage == 'full':
            prefix, env = pipeline_command(self.sources[version][1], args.python)
            cmd = [*prefix, '-i', str(inputs), '-o', str(out), '-c', str(args.cores),
                   '-t', *map(str, args.threshold), '-m', args.mode,
                   '-u', args.usalign, '-b', args.beem]
            row = run_measured(cmd, self.python_env(env), log, args)
            models = out / 'models'
        elif version == 'usalign_dir':
            models = self.prepared_models(size, inputs)
            chain_list = out / 'chain_list'
            chain_list.write_text(''.join(f'{p.name}\n' for p in sorted(models.glob('*.pdb'))),
                                  encoding='utf-8')
            table = out / 'usalign_dir_output.txt'
            cmd = [args.usalign, '-dir', f'{models}/', str(chain_list), *shlex.split(args.raw_opts)]
            row = run_measured(cmd, {}, log, args, stdout_path=table)
            row['align_s'] = row['wall_s']
            if row['status'] == 'ok' and '-outfmt 2' in ' '.join(shlex.split(args.raw_opts)):
                subprocess.run([args.python, str(DRIVER), 'rawtable', '--table', str(table),
                                '--models', str(models), '--out', str(out)], check=True)
        else:
            models = self.prepared_models(size, inputs)
            cmd = [args.python, '-u', str(DRIVER), 'align', '--src', str(self.sources[version][1]),
                   '--models', str(models), '--out', str(out), '--usalign', args.usalign,
                   '--cores', str(args.cores)]
            row = run_measured(cmd, self.python_env({}), log, args)

        matrix = out / 'results' / 'us-align_score_matrix.csv'
        if row['status'] == 'ok' and version != 'usalign_dir' and not matrix.exists():
            row['status'] = 'no_matrix'
        row['chains'], row['residues'] = chain_stats(models)
        cache = out / 'results' / 'usalign_cache'
        if cache.exists():
            cache_mb, row['cache_files'] = dir_usage(cache)
            row['cache_mb'] = round(cache_mb, 1)
        return row


# Results

def load_done(results_csv: Path):
    '''Rows already recorded, keyed by (stage, version, size, repeat).'''
    if not results_csv.exists():
        return {}
    with open(results_csv, newline='', encoding='utf-8') as f:
        return {(r.get('stage') or 'full', r['version'], int(r['size']), int(r['repeat'])): r
                for r in csv.DictReader(f)}


def append_row(results_csv: Path, row):
    '''Append one result row, writing the header for a new file.'''
    new = not results_csv.exists()
    with open(results_csv, 'a', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, extrasaction='ignore')
        if new:
            writer.writeheader()
        writer.writerow(row)


def compare_matrices(workdir: Path, stages, sizes):
    '''|TM difference| of every contender against the new version, per stage and size.'''
    try:
        import numpy as np  # pylint: disable=import-outside-toplevel
        import pandas as pd  # pylint: disable=import-outside-toplevel
    except ImportError:
        return
    rows = []
    for stage in stages:
        for size in sizes:
            ref_path = workdir / 'matrices' / f'{stage}_new_n{size:05d}.csv'
            if not ref_path.exists():
                continue
            ref = pd.read_csv(ref_path, index_col=0)
            for version in STAGE_CONTENDERS[stage]:
                path = workdir / 'matrices' / f'{stage}_{version}_n{size:05d}.csv'
                if version == 'new' or not path.exists():
                    continue
                other = pd.read_csv(path, index_col=0)
                common = ref.index.intersection(other.index)
                diff = (ref.loc[common, common] - other.loc[common, common]).abs().to_numpy()
                off = diff[~np.eye(len(common), dtype=bool)]
                off = off[~np.isnan(off)]
                rows.append({'stage': stage, 'size': size, 'versus_new': version,
                             'chains_compared': len(common),
                             'only_in_new': len(ref.index.difference(other.index)),
                             'only_in_other': len(other.index.difference(ref.index)),
                             'max_abs_diff': round(float(off.max()) if off.size else 0.0, 5),
                             'mean_abs_diff': round(float(off.mean()) if off.size else 0.0, 6),
                             'pairs_diff_gt_0.01': int((off > 0.01).sum() // 2)})
    if rows:
        out = workdir / 'matrix_agreement.csv'
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f'\nTM-score agreement with the new version ({out}):')
        print(pd.DataFrame(rows).to_string(index=False))


def summarise(workdir: Path, results_csv: Path):
    '''Print per-stage tables and save plots (needs pandas; plots need matplotlib).'''
    try:
        import pandas as pd  # pylint: disable=import-outside-toplevel
    except ImportError:
        print(f'Results in {results_csv}')
        return
    df = pd.read_csv(results_csv)
    if 'stage' not in df:
        df['stage'] = 'full'
    failed = df[df.status != 'ok'][['stage', 'version', 'size', 'repeat', 'status']]
    for stage, part in df[df.status == 'ok'].groupby('stage'):
        med = part.groupby(['version', 'size']).median(numeric_only=True)
        versions = [v for v in STAGE_CONTENDERS[stage] if v in med.index.get_level_values(0)]
        table = pd.DataFrame({'chains': med['chains'].groupby('size').max()})
        for v in versions:
            table[f'{v}_wall_s'] = med.loc[v, 'wall_s']
            table[f'{v}_cpu_s'] = med.loc[v, 'cpu_s']
            table[f'{v}_peak_main_mb'] = med.loc[v, 'peak_main_rss_mb']
            table[f'{v}_peak_tree_mb'] = med.loc[v, 'peak_tree_rss_mb']
        if 'new' in versions:
            for v in versions:
                if v != 'new':
                    table[f'speedup_vs_{v}'] = (table[f'{v}_wall_s'] / table['new_wall_s']).round(2)
                    table[f'cpu_ratio_{v}_to_new'] = (table[f'{v}_cpu_s'] / table['new_cpu_s']).round(2)
        out = workdir / f'summary_{stage}.csv'
        table.to_csv(out)
        print(f'\n[{stage}] median per size ({out}):')
        print(table.to_string())
        plot_stage(workdir, stage, med, versions, table)
    if not failed.empty:
        print('\nRuns that did not finish:')
        print(failed.to_string(index=False))


def plot_stage(workdir, stage, med, versions, table):
    '''Wall time, speed-up of the new version and peak memory for one stage.'''
    try:
        import matplotlib  # pylint: disable=import-outside-toplevel
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt  # pylint: disable=import-outside-toplevel
        from matplotlib.ticker import LogLocator, NullFormatter, ScalarFormatter  # pylint: disable=import-outside-toplevel
    except ImportError:
        return
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for v in versions:
        part = med.loc[v]
        colour, label = COLOURS[v], LABELS[v]
        axes[0].plot(part.index, part.wall_s, 'o-', color=colour, label=label)
        if v != 'new' and f'speedup_vs_{v}' in table:
            axes[1].plot(table.index, table[f'speedup_vs_{v}'], 'o-', color=colour,
                         label=f'vs {label}')
        axes[2].plot(part.index, part.peak_main_rss_mb, 'o-', color=colour, label=f'{label}, main')
        axes[2].plot(part.index, part.peak_tree_rss_mb, 's--', color=colour, label=f'{label}, all')
    axes[1].axhline(1, color='grey', lw=0.8)
    titles = ('Wall time', 'Speed-up of new version (other / new)', 'Peak memory')
    for ax, title, ylabel in zip(axes, titles, ('seconds', 'x', 'MB (summed RSS)')):
        ax.set_xscale('log')
        ax.set_xticks(list(table.index))
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis='x', labelrotation=45)
        ax.set_xlabel('structures in subset')
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    axes[2].legend(fontsize=7, ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.28))
    for ax in (axes[0], axes[2]):
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(LogLocator(subs=(1, 2, 5)))
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
    fig.suptitle({'full': 'Full pipeline', 'align': 'Alignment step only'}[stage])
    fig.tight_layout()
    path = workdir / f'benchmark_{stage}.png'
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'Plot: {path}')


def main():
    '''Run the benchmark.'''
    args = parse_args()
    workdir = args.workdir.resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    for tool in (args.usalign, args.beem, 'TreeCluster.py'):
        if shutil.which(tool) is None:
            sys.exit(f'{tool} not found on PATH')
    args.usalign, args.beem = shutil.which(args.usalign), shutil.which(args.beem)
    for folder in ('logs', 'matrices', 'runs'):
        workdir.joinpath(folder).mkdir(exist_ok=True)

    bench = Bench(args, workdir)
    subsets = make_subsets(args.structures, args.sizes, args.seed, workdir)
    results_csv = workdir / 'results.csv'
    done = load_done(results_csv)
    stopped = {(k[0], k[1]) for k, r in done.items() if r['status'] != 'ok'}

    for size, inputs in subsets.items():
        for repeat in range(args.repeats):
            for stage in args.stages:
                contenders = [v for v in STAGE_CONTENDERS[stage] if v in args.contenders]
                if not contenders:
                    continue
                shift = (size + repeat) % len(contenders)  # rotate run order
                for version in contenders[shift:] + contenders[:shift]:
                    key = (stage, version, size, repeat)
                    if key in done:
                        continue
                    base = {'stage': stage, 'version': version, 'ref': bench.refs[version],
                            'size': size, 'repeat': repeat}
                    if (stage, version) in stopped:
                        append_row(results_csv, {**base, 'status': 'skipped'})
                        continue
                    tag = f'{stage}_{version}_n{size:05d}_r{repeat}'
                    print(f'[{time.strftime("%H:%M:%S")}] {stage}: {LABELS[version]} n={size} '
                          f'repeat {repeat + 1}/{args.repeats} ...', flush=True)
                    out = workdir / 'runs' / tag
                    row = {**base, **bench.run(stage, version, size, inputs, out,
                                               workdir / 'logs' / f'{tag}.log')}
                    matrix = out / 'results' / 'us-align_score_matrix.csv'
                    if repeat == 0 and matrix.exists():
                        shutil.copy(matrix, workdir / 'matrices' / f'{stage}_{version}_n{size:05d}.csv')
                    if not args.keep_outputs:
                        shutil.rmtree(out, ignore_errors=True)
                    append_row(results_csv, row)
                    print(f'    {row["status"]}: {row["wall_s"]} s wall, {row["cpu_s"]} s CPU, '
                          f'{row["chains"]} chains, peak {row["peak_tree_rss_mb"]:.0f} MB', flush=True)
                    if row['status'] != 'ok':
                        print(f'    {LABELS[version]} stopped for the {stage} stage; see {row["log"]}')
                        stopped.add((stage, version))

    compare_matrices(workdir, args.stages, subsets)
    summarise(workdir, results_csv)


if __name__ == '__main__':
    main()
