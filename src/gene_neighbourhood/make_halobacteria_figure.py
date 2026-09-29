"""Single-panel summary of the gene neighbourhood analysis in Halobacteria.

For a set of anchor genes (one per genome), neighbour families within ±WINDOW are ranked by the
fraction of genomes in which they are found near the anchor. Cinquedea is compared with:
    - control gene families: random families present in a single copy in at least half of the
      Halobacteria genomes (what the neighbourhood of a conserved gene looks like);
    - random genes: one random gene per genome (what chance alone gives).

Run from the repo root:
    python src/gene_neighbourhood/make_halobacteria_figure.py [window in kb, default 10]
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

sys.path.insert(0, 'src/gene_neighbourhood')
from neighbourhood import (SEED, W, cinquedea_anchors, family_anchors, load_genes, load_homologs,
                           load_taxonomy, neighbourhoods)

font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if not f.fname.endswith('.ttc')]
for f in glob.glob(os.path.expanduser('~/.cache/fonts/Helvetica*.ttf')):
    font_manager.fontManager.addfont(f)
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 8, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False})

CLADE = 'Halobacteria'
COLOR = '#f58231'
N_RANKS = 20
N_CONTROLS = 50
N_RANDOM_DRAWS = 50


def neighbour_sharing(genes, anchor_ids, window):
    """Fraction of genomes in which each neighbour family is found near the anchor, most shared first."""
    nb, meta = neighbourhoods(genes, anchor_ids, window)
    p = nb.drop_duplicates(['anchor', 'family']).family.value_counts() / len(meta)
    # Ties broken by family identifier so that the ranking is reproducible
    return p.loc[sorted(p.index, key=lambda f: (-p[f], f))]


def sharing_curve(genes, anchor_ids, window):
    """Sharing values of the top N_RANKS neighbour families (padded with zeros)."""
    p = neighbour_sharing(genes, anchor_ids, window).values
    return np.pad(p[:N_RANKS], (0, max(0, N_RANKS - len(p))))


def compute_curves(genes, homologs, taxonomy, window):
    """Sharing curves of cinquedea, control families and random genes in Halobacteria genomes."""
    genomes = set(taxonomy.index[taxonomy.clade == CLADE])
    g = genes[genes.genome.isin(genomes)]
    cinquedea_families = set(genes.set_index('gene_id').loc[homologs.gene_id, 'family'])

    anchors = cinquedea_anchors(homologs)
    anchors = anchors[anchors.genome.isin(genomes)]
    cinquedea_neighbours = neighbour_sharing(g, anchors.gene_id, window)
    cinquedea = np.pad(cinquedea_neighbours.values[:N_RANKS], (0, max(0, N_RANKS - len(cinquedea_neighbours))))

    copies = g.dropna(subset=['family']).groupby(['family', 'genome']).size()
    single = copies[copies == 1].reset_index().groupby('family').genome.nunique()
    candidates = sorted(set(single.index[single >= 0.5 * len(genomes)]) - cinquedea_families)
    rng = np.random.default_rng(SEED)
    controls = rng.choice(candidates, N_CONTROLS, replace=False)
    control_curves = np.array([sharing_curve(g, family_anchors(g, f).gene_id, window) for f in controls])

    random_curves = np.array([
        sharing_curve(g, g.groupby('genome').sample(1, random_state=i).gene_id, window)
        for i in range(N_RANDOM_DRAWS)
    ])
    return {'cinquedea': cinquedea, 'cinquedea_neighbours': cinquedea_neighbours, 'controls': control_curves,
            'random': random_curves, 'n_genomes': len(anchors)}


def plot_curves(ax, curves, window_kb):
    x = np.arange(1, N_RANKS + 1)
    for key, color, label, ls in [
        ('controls', '#555555', f'Control gene families (n={N_CONTROLS})', '-'),
        ('random', '#aaaaaa', 'Random genes', '--'),
    ]:
        c = curves[key]
        ax.fill_between(x, np.quantile(c, 0.25, 0), np.quantile(c, 0.75, 0), color=color, alpha=0.2, lw=0)
        ax.plot(x, np.median(c, 0), color=color, ls=ls, lw=1.2, label=label)
    ax.plot(x, curves['cinquedea'], color=COLOR, marker='o', ms=3, lw=1.5, label='Cinquedea')

    ax.set_xlim(0.5, N_RANKS + 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([1, 5, 10, 15, 20])
    ax.set_xlabel('Neighbour family (ranked by conservation)')
    ax.set_ylabel(f'Fraction of genomes with the neighbour\nwithin ±{window_kb} kb')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], frameon=False, fontsize=7, loc='upper right')


def plot_neighbour_bars(ax, curves, window_kb, names, colors, n=10):
    """Cinquedea's most shared neighbour families (bars) against controls and random genes at each rank."""
    top = curves['cinquedea_neighbours'].iloc[:n]
    y = np.arange(len(top))
    ax.barh(y, top.values, color=[colors.get(f, '#cccccc') for f in top.index], edgecolor='black', lw=0.4,
            height=0.7, label='Cinquedea neighbours', zorder=2)

    c = curves['controls'][:, :n]
    q1, q2, q3 = np.quantile(c, 0.25, 0), np.median(c, 0), np.quantile(c, 0.75, 0)
    ax.errorbar(q2, y, xerr=[q2 - q1, q3 - q2], fmt='D', color='#333333', ms=3.5, lw=0.8, capsize=2,
                label=f'Control gene families, same rank (median, IQR; n={N_CONTROLS})', zorder=3)
    r = np.median(curves['random'][:, :n], 0)
    ax.plot(r, y, color='#999999', ls='--', lw=1, label='Random genes, same rank (median)', zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels([f'{names[f]}  ({v:.0%})' for f, v in top.items()], fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
    ax.set_xlabel(f'Genomes with the family within ±{window_kb} kb of the anchor gene')
    # Neutral proxy for the bars (their colours follow the gene families of panels A and B)
    from matplotlib.patches import Patch
    handles, labels = ax.get_legend_handles_labels()
    handles = [Patch(facecolor='#cccccc', edgecolor='black', lw=0.4) if l == 'Cinquedea neighbours' else h
               for h, l in zip(handles, labels)]
    order = [labels.index('Cinquedea neighbours')] + [i for i, l in enumerate(labels) if l.startswith('Control')] + \
            [i for i, l in enumerate(labels) if l.startswith('Random')]
    ax.legend([handles[i] for i in order], [labels[i] for i in order], frameon=False, fontsize=6.5,
              loc='upper left', bbox_to_anchor=(-0.02, -0.22))


def save_curves(curves, window_kb):
    c, r = curves['controls'], curves['random']
    pd.DataFrame({
        'rank': np.arange(1, N_RANKS + 1),
        'cinquedea': curves['cinquedea'],
        'controls_median': np.median(c, 0),
        'controls_q25': np.quantile(c, 0.25, 0),
        'controls_q75': np.quantile(c, 0.75, 0),
        'random_median': np.median(r, 0),
        'random_q25': np.quantile(r, 0.25, 0),
        'random_q75': np.quantile(r, 0.75, 0),
    }).to_csv(W / f'window_{window_kb}kb' / 'halobacteria_sharing_curves.tsv', sep='\t', index=False)


if __name__ == '__main__':
    window_kb = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    genes = load_genes()
    taxonomy = load_taxonomy(sorted(genes.genome.unique()))
    curves = compute_curves(genes, load_homologs(), taxonomy, window_kb * 1000)
    save_curves(curves, window_kb)

    fig, ax = plt.subplots(figsize=(3.5, 2.8))
    plot_curves(ax, curves, window_kb)
    ax.set_title(f'{CLADE} (n={curves["n_genomes"]} genomes)', loc='left', fontsize=8)

    fig_dir = W / 'figures'
    fig.savefig(fig_dir / f'fig_neighbourhood_halobacteria_{window_kb}kb.pdf', bbox_inches='tight')
    fig.savefig(fig_dir / f'fig_neighbourhood_halobacteria_{window_kb}kb.png', bbox_inches='tight', dpi=300)

    top = curves['controls'][:, 0]
    print(f'±{window_kb} kb: cinquedea top neighbour {curves["cinquedea"][0]:.2f}; controls median {np.median(top):.2f} '
          f'(cinquedea below {np.mean(top > curves["cinquedea"][0]):.0%} of controls); '
          f'random {np.median(curves["random"][:, 0]):.2f}')
