'''Panel-size targeting and panel assembly.

A panel is one structure per cluster (the cluster representative, or the
structure itself for a singleton) plus a few low-coverage additions: the
structures least well represented by the cluster panel, added farthest-first.

Choosing a threshold for a requested panel size reproduces TreeCluster's
max_clade rule on the written Newick tree. Every clade has a maximum pairwise
leaf distance (its diameter); the clusters at threshold t are the largest clades
with a diameter of at most t. Diameters never decrease towards the root, so the
cluster count falls by (children - 1) each time t reaches an internal node's
diameter. Sweeping those diameters gives every cluster count the tree can
produce, without running TreeCluster once per candidate threshold.
'''

import csv
from dataclasses import dataclass, field
from math import floor, inf

import numpy as np
from Bio import Phylo

from .common import extract_sequence

DEFAULT_LOW_COVERAGE_FRACTION = 0.05
MAX_LOW_COVERAGE_FRACTION = 0.5


def parse_panel_size(text):
    '''Parse "N" or "MIN-MAX" into an inclusive (min, max) tuple of positive integers.'''
    parts = str(text).strip().split('-')
    try:
        if len(parts) == 1:
            low = high = int(parts[0])
        elif len(parts) == 2:
            low, high = int(parts[0]), int(parts[1])
        else:
            raise ValueError
    except ValueError:
        raise ValueError(
            f'Invalid panel size {text!r}: use a whole number (e.g. 50) '
            f'or an inclusive range (e.g. 40-60).'
        ) from None
    if low < 1 or high < low:
        raise ValueError(
            f'Invalid panel size {text!r}: sizes must be at least 1, '
            f'and a range must be written MIN-MAX.'
        )
    return low, high


def format_threshold(threshold):
    '''Two decimals where exact (0.20 -> "0.20"), otherwise as many as needed.'''
    threshold = float(threshold)
    for decimals in range(2, 13):
        text = f'{threshold:.{decimals}f}'
        if float(text) == threshold:
            return text
    return repr(threshold)


def low_coverage_count(n_clusters, fraction, pool_size):
    '''Additions so that they make up `fraction` of the final panel (rounded half up).'''
    if fraction <= 0:
        return 0
    wanted = floor(fraction * n_clusters / (1 - fraction) + 0.5)
    return max(0, min(wanted, pool_size - n_clusters))


def _clade_diameters(tree_path):
    '''Leaf count and the diameter of every internal node with its number of children.

    Mirrors TreeCluster.min_clusters_threshold_max_clade, including its order of
    floating-point additions, so threshold comparisons agree with TreeCluster.
    '''
    tree = Phylo.read(str(tree_path), 'newick')
    leaf_dist, max_pair, internal = {}, {}, []
    for node in tree.find_clades(order='postorder'):
        if node.is_terminal():
            leaf_dist[node] = 0.0
            max_pair[node] = 0.0
            continue
        best = second = -inf
        for child in node.clades:
            dist = leaf_dist[child] + (child.branch_length or 0.0)
            if dist > best:
                second, best = best, dist
            elif dist > second:
                second = dist
        leaf_dist[node] = best
        max_pair[node] = max([max_pair[c] for c in node.clades] + [best + second])
        internal.append((max_pair[node], len(node.clades)))
    return len(tree.get_terminals()), internal


@dataclass
class CutOption:
    '''One achievable clustering of the tree.'''
    n_clusters: int
    threshold_min: float      # smallest threshold that gives n_clusters
    threshold_max: float      # thresholds below this value give n_clusters
    n_low_coverage: int
    requested: list = field(default_factory=list)

    @property
    def panel_size(self):
        '''Clusters plus low-coverage additions.'''
        return self.n_clusters + self.n_low_coverage

    @property
    def threshold(self):
        '''The shortest decimal in [threshold_min, threshold_max), as near its middle
        as that precision allows, so that it is not on the edge of a merge.'''
        low, high = self.threshold_min, self.threshold_max
        middle = low + 0.01 if high == inf else (low + high) / 2
        for decimals in range(2, 13):
            scale = 10 ** decimals
            first = int(np.ceil(low * scale - 1e-9))
            last = first if high == inf else int(np.ceil(high * scale - 1e-9)) - 1
            if first > last:
                continue
            value = round(min(max(round(middle * scale), first), last) / scale, decimals)
            if low <= value < high:
                return value
        return low


def cut_options(tree_path, low_coverage_fraction=DEFAULT_LOW_COVERAGE_FRACTION):
    '''Every distinct cluster count the tree can give, smallest threshold first.'''
    n_leaves, internal = _clade_diameters(tree_path)
    drops = {}
    for diameter, n_children in internal:
        drops[diameter] = drops.get(diameter, 0) + n_children - 1
    options, n_clusters, start = [], n_leaves, 0.0
    for diameter in sorted(drops):
        if diameter > start:
            options.append((n_clusters, start, diameter))
        n_clusters -= drops[diameter]
        start = max(start, diameter)
    options.append((n_clusters, start, inf))
    return [
        CutOption(k, low, high, low_coverage_count(k, low_coverage_fraction, n_leaves))
        for k, low, high in options if k >= 1
    ]


def choose_cuts(options, panel_sizes):
    '''Pick the cuts for each requested (min, max) panel size.

    Every cut whose panel size lies in the range is kept. If none does, the cut
    closest to the range is used, preferring the smaller panel on a tie.
    Where several cuts give the same panel size, the one with most clusters is kept.
    '''
    chosen = {}
    for low, high in panel_sizes:
        label = str(low) if low == high else f'{low}-{high}'
        hits = [o for o in options if low <= o.panel_size <= high]
        if not hits:
            def gap(option):
                size = option.panel_size
                return (low - size if size < low else size - high, size)
            hits = [min(options, key=gap)]
        # Near the top of the tree, the pool-size cap can give several cuts the same
        # panel size; keep the one with the most clusters (fewest additions).
        by_size = {}
        for option in hits:
            best = by_size.get(option.panel_size)
            if best is None or option.n_clusters > best.n_clusters:
                by_size[option.panel_size] = option
        for option in by_size.values():
            option.requested.append(label)
            chosen[option.n_clusters] = option
    return sorted(chosen.values(), key=lambda o: o.n_clusters)


def write_cut_options(options, path):
    '''Write every achievable panel size and the thresholds that produce it.'''
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow([
            'panel_size', 'clusters', 'low_coverage_additions',
            'threshold_from', 'threshold_below', 'threshold_used', 'requested',
        ])
        for option in sorted(options, key=lambda o: o.n_clusters, reverse=True):
            used = format_threshold(option.threshold) if option.requested else ''
            writer.writerow([
                option.panel_size, option.n_clusters, option.n_low_coverage,
                f'{option.threshold_min:.6g}',
                '' if option.threshold_max == inf else f'{option.threshold_max:.6g}',
                used, ';'.join(option.requested),
            ])


def add_low_coverage(tmmatrix, panel, n_additions):
    '''Farthest-first additions: repeatedly add the structure whose best TM-score
    to the current panel is lowest. Ties go to the first name alphabetically.
    Returns [(name, best TM-score to the panel before it was added)].'''
    names = sorted(tmmatrix.index)
    scores = tmmatrix.loc[names, names].to_numpy(dtype=float)
    index = {name: i for i, name in enumerate(names)}
    members = [index[name] for name in panel]
    coverage = scores[:, members].max(axis=1) if members else np.full(len(names), -inf)
    available = np.ones(len(names), dtype=bool)
    available[members] = False
    added = []
    while len(added) < n_additions and available.any():
        pick = int(np.argmin(np.where(available, coverage, inf)))
        added.append((names[pick], float(coverage[pick])))
        available[pick] = False
        coverage = np.maximum(coverage, scores[:, pick])
    return added


def _medoid(tmmatrix, members):
    sub = tmmatrix.loc[list(members), list(members)]
    return sub.mean(axis=0).idxmax()


# pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
def build_panel(threshold_str, clusters, representatives, tmmatrix, low_coverage_fraction,
                models_by_name, out_dir, requested=None):
    '''Assemble, score and write one threshold's panel. Returns a summary row.'''
    rows = []
    for cluster_num in sorted(clusters):
        members = clusters[cluster_num]
        if cluster_num == -1:
            rows.extend({'structure': name, 'role': 'singleton', 'cluster': -1,
                         'cluster_size': 1} for name in sorted(members))
            continue
        entry = representatives.get(cluster_num)
        if entry is not None:
            name, role = entry['record_name'], 'cluster_representative'
        else:
            name, role = _medoid(tmmatrix, members), 'cluster_medoid'
            print(f'Warning: no consensus representative for cluster {cluster_num} at '
                  f'threshold {threshold_str}; using its medoid ({name}) in the panel.')
        rows.append({'structure': name, 'role': role, 'cluster': cluster_num,
                     'cluster_size': len(members)})

    n_clusters = len(rows)
    n_low = low_coverage_count(n_clusters, low_coverage_fraction, len(tmmatrix))
    cluster_of = {name: num for num, members in clusters.items() for name in members}
    for name, before in add_low_coverage(tmmatrix, [r['structure'] for r in rows], n_low):
        cluster = cluster_of[name]  # never a singleton: those are already in the panel
        rows.append({'structure': name, 'role': 'low_coverage', 'cluster': cluster,
                     'cluster_size': len(clusters[cluster]), 'tm_to_panel_before': before})

    panel = [r['structure'] for r in rows]
    best = tmmatrix.loc[:, panel].max(axis=1)
    nearest = tmmatrix.loc[:, panel].idxmax(axis=1)
    represented = nearest.value_counts()
    for row in rows:
        row['structures_represented'] = int(represented.get(row['structure'], 0))

    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f'threshold_{threshold_str}_panel'
    columns = ['structure', 'role', 'cluster', 'cluster_size', 'structures_represented',
               'tm_to_panel_before']
    with open(out_dir.joinpath(f'{stem}.csv'), 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, restval='')
        writer.writeheader()
        writer.writerows(rows)
    with open(out_dir.joinpath(f'{stem}.fasta'), 'w', encoding='utf-8') as handle:
        for row in rows:
            sequence = extract_sequence(models_by_name[row['structure']])
            handle.write(f'>{row["structure"]} role={row["role"]} cluster={row["cluster"]}\n')
            for i in range(0, len(sequence), 80):
                handle.write(f'{sequence[i:i + 80]}\n')

    return {
        'threshold': threshold_str,
        'requested_panel_size': ';'.join(requested or []),
        'panel_size': len(rows),
        'clusters': sum(r['role'] in ('cluster_representative', 'cluster_medoid') for r in rows),
        'singletons': sum(r['role'] == 'singleton' for r in rows),
        'low_coverage_additions': len(rows) - n_clusters,
        'mean_best_tm': float(best.mean()),
        'p10_best_tm': float(np.quantile(best, 0.1)),
        'min_best_tm': float(best.min()),
        'panel_csv': f'clusters/{threshold_str}/{stem}.csv',
    }


def write_panel_summary(rows, path):
    '''One row per threshold: panel composition and how well it covers every structure.'''
    if not rows:
        return
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v)
                             for k, v in row.items()})
