"""Figures for the gene neighbourhood analysis of cinquedea homologs.

Run from the repo root (after run_analysis.py and predict_pfam.sh):
    python src/gene_neighbourhood/make_figures.py [window in kb, default 10]

Outputs (data/gene_neighbourhood/figures/):
    fig_neighbourhood_<N>kb.pdf        (A) neighbours shared across homologs, (B) conservation vs
                                       phylogenetic distance, cinquedea vs control families
    fig_neighbourhood_top10_<N>kb.pdf            gene maps around the top 10 homologs, one per species
                                                 (addendum to Figure 3C)
    fig_neighbourhood_top10_per_genus_<N>kb.pdf  same, one per genus
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dna_features_viewer import GraphicFeature, GraphicRecord
from matplotlib import font_manager
from matplotlib.patches import Patch

sys.path.insert(0, 'src/gene_neighbourhood')
sys.path.insert(0, 'src/utils')
from neighbourhood import QUERY_GENOME, W, load_family_pfam, load_genes, load_homologs
from palette import palette_20

# Helvetica, to match the other figures of the paper (see src/foldseek_2090/make_figures.py).
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if not f.fname.endswith('.ttc')]
for f in glob.glob(os.path.expanduser('~/.cache/fonts/Helvetica*.ttf')):
    font_manager.fontManager.addfont(f)
plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 8, 'pdf.fonttype': 42,
                     'axes.spines.top': False, 'axes.spines.right': False})

CLADE_COLORS = {
    'Halobacteria': '#f58231',
    'Nitrososphaeria': '#9a6324',
    'Other archaea': '#ffd8b1',
    'Bacteria': '#4363d8',
}
CLADE_ORDER = ['Halobacteria', 'Nitrososphaeria', 'Other archaea', 'Bacteria']
CINQUEDEA_COLOR = '#3cb44b'
N_TOP_FAMILIES = 6
FIG = W / 'figures'


def family_label(family, pfam):
    arch = pfam.get(family)
    return arch.split('+')[0] if isinstance(arch, str) else 'no Pfam'


def relative_positions(nb, genes):
    """Gene position relative to the anchor (in genes, oriented by the anchor's strand)."""
    g = genes.set_index('gene_id')
    a = g.loc[nb.anchor]
    return np.where(a.strand.values == '+', nb.gene_index.values - a.gene_index.values,
                    a.gene_index.values - nb.gene_index.values)


def panel_heatmap(ax, nb, anchors, enrichment, pfam):
    top = (enrichment[enrichment.clade.isin(['Halobacteria', 'Nitrososphaeria']) & (enrichment.q < 0.05)]
           .groupby('clade').head(N_TOP_FAMILIES))
    families = list(dict.fromkeys(top.family))
    pos = nb.groupby('family').rel.agg(lambda s: s.mode().iloc[0])

    a = anchors.assign(clade=pd.Categorical(anchors.clade, CLADE_ORDER)).sort_values(
        ['clade', 'order', 'family', 'genus', 'species']).reset_index(drop=True)
    present = nb[nb.family.isin(families)].groupby(['anchor', 'family']).size().unstack(fill_value=0) > 0
    present = present.reindex(index=a.anchor, columns=families, fill_value=False)

    img = np.ones((len(a), len(families), 3))
    for i, clade in enumerate(a.clade):
        rgb = plt.matplotlib.colors.to_rgb(CLADE_COLORS[clade])
        img[i, present.iloc[i].values] = rgb
    ax.imshow(img, aspect='auto', interpolation='nearest')
    ax.set_xticks(range(len(families)))
    ax.set_xticklabels([f'{family_label(f, pfam)} ({pos[f]:+d})' for f in families], rotation=60, ha='right', fontsize=6.5)
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    # Clade separators, side bar and labels
    bounds = a.groupby('clade', observed=True).apply(lambda d: (d.index.min(), d.index.max()))
    for clade, (lo, hi) in bounds.items():
        ax.axhline(hi + 0.5, color='black', lw=0.4)
        ax.add_patch(plt.Rectangle((-1.1, lo - 0.5), 0.5, hi - lo + 1, color=CLADE_COLORS[clade],
                                   clip_on=False, transform=ax.transData))
        if hi - lo >= 10:
            ax.text(-1.3, (lo + hi) / 2, f'{clade}\n(n={hi - lo + 1})', ha='right', va='center', fontsize=6.5)
    small = [f'{c} (n={hi - lo + 1})' for c, (lo, hi) in bounds.items() if hi - lo < 10]
    ax.text(-1.3, len(a) - 0.5, '\n'.join(small), ha='right', va='top', fontsize=6)
    ax.set_xlim(-1.1, len(families) - 0.5)


def panel_conservation(ax, pairs, clade, bins):
    d = pairs[(pairs.clade_a == clade) & (pairs.clade_b == clade)].copy()
    d['bin'] = pd.cut(d.distance, bins)
    x = np.arange(len(bins) - 1)
    med = d.groupby(['anchor_family', 'bin']).jaccard.median().unstack()
    controls = med.drop(index='cinquedea')
    for _, row in controls.iterrows():
        ax.plot(x, row.values, color='#bbbbbb', lw=0.6, zorder=1)
    ax.plot(x, controls.median().values, color='black', lw=1.2, ls='--', label='Control families (median)', zorder=2)

    c = d[d.anchor_family == 'cinquedea'].groupby('bin').jaccard
    q1, q2, q3 = c.quantile(0.25).values, c.median().values, c.quantile(0.75).values
    ax.errorbar(x, q2, yerr=[q2 - q1, q3 - q2], color=CLADE_COLORS[clade], marker='o', ms=4, lw=1.4,
                capsize=2, label='Cinquedea (median, IQR)', zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{max(b.left, 0):g}–{b.right:g}' for b in med.columns], rotation=45, ha='right', fontsize=6.5)
    ax.set_ylim(0, 1)
    ax.set_title(clade, fontsize=8, color=CLADE_COLORS[clade])
    ax.set_xlabel('Phylogenetic distance (GTDB r214 tree)')


def figure_main(window_kb, genes, pfam):
    d = W / f'window_{window_kb}kb'
    nb = pd.read_csv(d / 'neighbourhood_genes.tsv.gz', sep='\t')
    nb['rel'] = relative_positions(nb, genes)
    anchors = pd.read_csv(d / 'anchors.tsv', sep='\t')
    enrichment = pd.read_csv(d / 'family_enrichment.tsv', sep='\t')
    pairs = pd.read_csv(d / 'conservation_pairs.tsv.gz', sep='\t')

    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.1, 1], height_ratios=[1, 1], wspace=0.35, hspace=0.9)
    ax = fig.add_subplot(gs[:, 0])
    panel_heatmap(ax, nb, anchors, enrichment, pfam)
    ax.set_title(f'A  Recurrent neighbours (±{window_kb} kb)', loc='left', fontweight='bold', x=-0.35)

    bins = [0, 0.025, 0.05, 0.1, 0.2, 0.4]
    ax1 = fig.add_subplot(gs[0, 1])
    panel_conservation(ax1, pairs, 'Halobacteria', bins)
    ax1.set_ylabel('Neighbourhood similarity\n(Jaccard)')
    ax1.set_title('B  Neighbourhood conservation', loc='left', fontweight='bold', x=-0.3, y=1.12)
    ax1.legend(frameon=False, fontsize=6.5, loc='upper right')
    ax2 = fig.add_subplot(gs[1, 1])
    panel_conservation(ax2, pairs, 'Nitrososphaeria', bins)
    ax2.set_ylabel('Neighbourhood similarity\n(Jaccard)')

    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / f'fig_neighbourhood_{window_kb}kb.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'fig_neighbourhood_{window_kb}kb.png', bbox_inches='tight', dpi=300)


def select_homologs(homologs, by, n):
    """Top homologs by search score, one per GTDB species or genus (cinquedea first)."""
    h = homologs.dropna(subset=['gene_id']).copy()
    h['search_bitscore'] = h.search_bitscore.fillna(np.inf)
    h['species'] = h.gtdb_species.fillna('Haloferax larsenii s5a-1')
    h['genus'] = h.species.str.split(' ').str[0]
    return h.sort_values('search_bitscore', ascending=False).drop_duplicates(by).head(n)


def figure_top(window_kb, genes, pfam, n=10, by='species', n_colors=12):
    top = select_homologs(load_homologs(), by, n)
    g = genes.set_index('gene_id')
    windows = []
    for _, h in top.iterrows():
        a = g.loc[h.gene_id]
        c = genes[(genes.contig == a.contig) & (genes.end >= a.start - window_kb * 1000) &
                  (genes.start <= a.end + window_kb * 1000)].sort_values('start')
        windows.append((h, a, c))

    # Colour the families shared by most of the displayed neighbourhoods
    anchor_families = set(top.gene_id.map(g.family))
    counts = pd.Series([f for _, _, c in windows for f in set(c.family.dropna()) - anchor_families]).value_counts()
    shared = list(counts[counts >= 2].index[:n_colors])
    palette = [c for c in palette_20 if c not in (CINQUEDEA_COLOR, '#808080', '#fffac8', '#ffe119')]
    colors = {f: palette[i] for i, f in enumerate(shared)}

    fig, axes = plt.subplots(len(windows), 1, figsize=(7.2, 0.55 * len(windows)))
    span = window_kb * 1000
    for ax, (h, a, c) in zip(axes, windows):
        flip = a.strand == '-'
        coords = []
        for _, x in c.iterrows():
            start, end = (a.end - x.end, a.end - x.start) if flip else (x.start - a.start, x.end - a.start)
            coords.append((start, end, (1 if x.strand == '+' else -1) * (-1 if flip else 1), x.gene_id, x.family))
        coords.sort()
        features = []
        for i, (start, end, strand, gene_id, family) in enumerate(coords):
            # Trim small overlaps between adjacent genes so that all genes are drawn on one line
            if i + 1 < len(coords):
                end = min(end, coords[i + 1][0] - 1)
            color = CINQUEDEA_COLOR if gene_id == h.gene_id else colors.get(family, '#e6e6e6')
            features.append(GraphicFeature(start=start, end=end, strand=strand, color=color,
                                           linewidth=0.4, thickness=9))
        length = a.end - a.start
        record = GraphicRecord(first_index=-span, sequence_length=2 * span + length, features=features)
        record.plot(ax=ax, with_ruler=False, draw_line=True)
        ax.set_xlim(-span, span + length)
        ax.set_ylim(-1, 1)
        label = h.species if by == 'species' else f'{h.species}'
        ax.text(-span, 0.95, label, fontsize=7, fontstyle='italic', va='bottom', transform=ax.get_xaxis_transform())
        ax.text(span + length, 0.95, 'this study' if h.genome == QUERY_GENOME else h.genome, fontsize=6,
                color='grey', ha='right', va='bottom', transform=ax.get_xaxis_transform())

    handles = [Patch(color=CINQUEDEA_COLOR, label='Cinquedea homolog')]
    unannotated = 0
    for f in shared:
        arch = pfam.get(f)
        if isinstance(arch, str):
            # First two Pfam domains of the representative (N- to C-terminal)
            label = ' + '.join(arch.split('+')[:2])
        else:
            unannotated += 1
            label = f'Unannotated family {unannotated}'
        handles.append(Patch(color=colors[f], label=label))
    handles.append(Patch(color='#e6e6e6', label='Other'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.02), ncol=3, frameon=False, fontsize=6.5)

    name = f'fig_neighbourhood_top{n}_{"per_genus_" if by == "genus" else ""}{window_kb}kb'
    fig.savefig(FIG / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'{name}.png', bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    window_kb = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    genes = load_genes()
    pfam = load_family_pfam()
    figure_main(window_kb, genes, pfam)
    figure_top(window_kb, genes, pfam, by='species')
    figure_top(window_kb, genes, pfam, by='genus')
