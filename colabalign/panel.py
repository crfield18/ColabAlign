'''Panel-size targeting and panel assembly.

A panel is the representative of every cluster with two or more members,
plus optional low-coverage additions, added farthest-first: each is the
structure least similar to the panel so far. Additions come from the
singletons (cluster -1) first; once every singleton is in the panel, any
remaining additions come from the other members of the clusters.

Choosing a threshold for a requested panel size reproduces TreeCluster's
max_clade rule on the written Newick tree. Every clade has a maximum pairwise
leaf distance (its diameter); the clusters at threshold t are the largest clades
with a diameter of at most t. Diameters never decrease towards the root, so as
t reaches an internal node's diameter, that node's children merge into one
cluster: its leaf children stop being singletons, and its internal children
(each already a cluster) become one. Sweeping those diameters gives every
clustering the tree can produce, without running TreeCluster for each one.
'''

import csv
from dataclasses import dataclass, field
from math import floor, inf

import numpy as np
import pandas as pd
from Bio import Phylo

from .common import extract_sequence

DEFAULT_LOW_COVERAGE_FRACTION = 0.0
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


def low_coverage_count(n_clusters, fraction, n_available):
    '''Additions so that they make up `fraction` of the final panel (rounded
    half up), limited by the structures not already in the panel.'''
    if fraction <= 0:
        return 0
    wanted = floor(fraction * n_clusters / (1 - fraction) + 0.5)
    return max(0, min(wanted, n_available))


def _clade_diameters(tree_path):
    '''Leaf count, and for every internal node its diameter, number of leaf
    children and number of internal children.

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
        leaves = sum(1 for c in node.clades if c.is_terminal())
        internal.append((max_pair[node], leaves, len(node.clades) - leaves))
    return len(tree.get_terminals()), internal


@dataclass
class CutOption:  # pylint: disable=too-many-instance-attributes
    '''One achievable clustering of the tree.'''
    n_clusters: int           # clusters with two or more members
    n_singletons: int         # cluster -1
    threshold_min: float      # smallest threshold that gives this clustering
    threshold_max: float      # thresholds below this value give it
    n_low_coverage: int       # low-coverage additions (singletons first)
    requested: list = field(default_factory=list)

    @property
    def n_groups(self):
        '''Clusters plus singletons: TreeCluster's group count.'''
        return self.n_clusters + self.n_singletons

    @property
    def panel_size(self):
        '''Cluster representatives plus low-coverage additions.'''
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
    '''Every clustering the tree can give, smallest threshold first, from the
    threshold with the most multi-member clusters upwards.

    Below that threshold, raising it pairs up singletons; above it, clusters
    merge. The same panel size can come from either side, but below the peak
    most structures are singletons, so those cuts are not offered.'''
    n_leaves, internal = _clade_diameters(tree_path)
    events = {}
    for diameter, leaf_children, internal_children in internal:
        clusters, singletons = events.get(diameter, (0, 0))
        events[diameter] = (clusters + 1 - internal_children, singletons - leaf_children)
    cuts, clusters, singletons, start = [], 0, n_leaves, 0.0
    for diameter in sorted(events):
        if diameter > start:
            cuts.append((clusters, singletons, start, diameter))
        clusters += events[diameter][0]
        singletons += events[diameter][1]
        start = max(start, diameter)
    cuts.append((clusters, singletons, start, inf))
    peak = max(range(len(cuts)), key=lambda i: (cuts[i][0], -i))
    return [
        CutOption(k, s1, low, high, low_coverage_count(k, low_coverage_fraction, n_leaves - k))
        for k, s1, low, high in cuts[peak:]
    ]


def choose_cuts(options, panel_sizes):
    '''Pick the cuts for each requested (min, max) panel size.

    Every cut whose panel size lies in the range is kept. If none does, the cut
    closest to the range is used, preferring the smaller panel on a tie. Where
    several cuts give the same panel size, the one with fewest singletons (the
    highest threshold) is kept, so the fewest structures are left unclustered.'''
    chosen = {}
    for low, high in panel_sizes:
        label = str(low) if low == high else f'{low}-{high}'
        hits = [o for o in options if low <= o.panel_size <= high]
        if not hits:
            def gap(option):
                size = option.panel_size
                return (low - size if size < low else size - high, size, option.n_singletons)
            hits = [min(options, key=gap)]
        by_size = {}
        for option in hits:
            best = by_size.get(option.panel_size)
            if best is None or option.n_singletons < best.n_singletons:
                by_size[option.panel_size] = option
        for option in by_size.values():
            option.requested.append(label)
            chosen[option.threshold_min] = option
    return sorted(chosen.values(), key=lambda o: o.threshold_min)


def write_cut_options(options, path):
    '''Write every achievable panel size and the thresholds that produce it.'''
    with open(path, 'w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow([
            'panel_size', 'clusters', 'low_coverage_additions', 'singletons',
            'threshold_from', 'threshold_below', 'threshold_used', 'requested',
        ])
        for option in options:
            used = format_threshold(option.threshold) if option.requested else ''
            writer.writerow([
                option.panel_size, option.n_clusters, option.n_low_coverage, option.n_singletons,
                f'{option.threshold_min:.6g}',
                '' if option.threshold_max == inf else f'{option.threshold_max:.6g}',
                used, ';'.join(option.requested),
            ])


def add_low_coverage(tmmatrix, panel, n_additions, candidates):
    '''Farthest-first additions from `candidates`: repeatedly add the candidate
    whose best TM-score to the current panel is lowest. Ties go to the first
    name alphabetically. With an empty panel, the first pick is the first name.
    Returns [(name, best TM-score to the panel before it was added)].'''
    names = sorted(tmmatrix.index)
    scores = tmmatrix.loc[names, names].to_numpy(dtype=float)
    index = {name: i for i, name in enumerate(names)}
    members = [index[name] for name in panel]
    coverage = scores[:, members].max(axis=1) if members else np.full(len(names), -inf)
    available = np.zeros(len(names), dtype=bool)
    available[[index[name] for name in candidates]] = True
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
    '''Assemble, score and write one threshold's panel. Returns a summary row.

    The panel is the representative of every cluster with two or more members,
    plus low-coverage additions: singletons (cluster -1) first, then, once all
    singletons are in, the least-covered members of the clusters.'''
    rows = []
    for cluster_num in sorted(c for c in clusters if c != -1):
        members = clusters[cluster_num]
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
    singletons = sorted(clusters.get(-1, []))
    n_low = low_coverage_count(n_clusters, low_coverage_fraction, len(tmmatrix) - n_clusters)
    cluster_of = {name: num for num, members in clusters.items() for name in members}
    # Singletons first; any slots left over go to the least-covered cluster members.
    for pool in (singletons, sorted(tmmatrix.index)):
        panel = [r['structure'] for r in rows]
        wanted = n_low - (len(rows) - n_clusters)
        for name, before in add_low_coverage(tmmatrix, panel, wanted, candidates=pool):
            cluster = cluster_of[name]
            rows.append({'structure': name, 'role': 'low_coverage', 'cluster': cluster,
                         'cluster_size': 1 if cluster == -1 else len(clusters[cluster]),
                         'tm_to_panel_before': before})

    panel = [r['structure'] for r in rows]
    if panel:
        best = tmmatrix.loc[:, panel].max(axis=1)
        represented = tmmatrix.loc[:, panel].idxmax(axis=1).value_counts()
    else:
        print(f'Warning: threshold {threshold_str} gives no clusters with two or more '
              f'members, so its panel is empty.')
        best, represented = pd.Series(np.nan, index=tmmatrix.index), pd.Series(dtype=int)
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
        'clusters': n_clusters,
        'low_coverage_additions': len(rows) - n_clusters,
        'singletons': len(singletons),
        'singletons_in_panel': sum(r['cluster'] == -1 for r in rows),
        'mean_best_tm': float(best.mean()),
        'p10_best_tm': float(best.quantile(0.1)),
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
