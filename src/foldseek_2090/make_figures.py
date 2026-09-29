"""Figures for the Foldseek analysis of hydrolase 2090 (cinquedea).

Run from the repo root:
    python src/foldseek_2090/make_figures.py

Figures are made twice: with all hits, and without eukaryotic hits (PCoA recomputed).
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd

sys.path.insert(0, 'src/foldseek_2090')
from foldseek_2090 import (DOMAIN_COLORS, PFAM_COLORS, PFAM_ORDER, QUERY, W,
                           annotate_hits, load_all_vs_all, load_pfam,
                           load_search_hits, pcoa, similarity_matrix)

# Helvetica, to match the other figures of the paper. Matplotlib cannot embed fonts from the
# macOS .ttc collection in PDFs, so individual faces were extracted to ~/.cache/fonts with
# fontTools (TTCollection('/System/Library/Fonts/Helvetica.ttc')).
font_manager.fontManager.ttflist = [f for f in font_manager.fontManager.ttflist if not f.fname.endswith('.ttc')]
for f in glob.glob(os.path.expanduser('~/.cache/fonts/Helvetica*.ttf')):
    font_manager.fontManager.addfont(f)

plt.rcParams.update({'font.family': 'Helvetica', 'font.size': 8, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})

DOMAIN_COLORS = {**DOMAIN_COLORS, 'Metagenome (unassigned)': '#008080'}
DOMAIN_ORDER = ['Eukaryota', 'Metagenome (unassigned)', 'Bacteria', 'Archaea']
FIG = f'{W}/figures'


def build_table(ava, exclude=(), suffix=''):
    hits = annotate_hits(load_search_hits(), load_pfam())
    hits['domain'] = hits.domain.fillna('Metagenome (unassigned)')
    hits = hits[~hits.domain.isin(exclude)]

    ids = sorted(set(ava['query']) & (set(hits.id) | {QUERY}))
    sim = similarity_matrix(ava, ids, 'mean_tm').clip(0, 1)
    coords, explained = pcoa(1 - sim, k=3)

    t = pd.DataFrame(coords, columns=['PC1', 'PC2', 'PC3'], index=ids)
    t['tm_to_2090'] = sim[ids.index(QUERY)]
    t = t.join(hits.set_index('id')[['taxname', 'domain', 'phylum', 'class', 'qtmscore', 'fident',
                                     'rank', 'pfam_domains', 'category']])
    t.loc[QUERY, ['taxname', 'domain', 'phylum', 'class', 'pfam_domains', 'category']] = [
        'Haloferax larsenii s5a-1', 'Archaea', 'Halobacteriota', 'Halobacteria', 'DUF900', 'DUF900']
    t.index.name = 'id'

    # Hits whose UniProt entry was retired have no AFDB v6 model, hence no PCoA coordinates.
    retired = hits[~hits.id.isin(ids)].set_index('id')
    retired = retired[['taxname', 'domain', 'phylum', 'class', 'qtmscore', 'fident', 'rank', 'pfam_domains', 'category']]
    table = pd.concat([t, retired])
    table['has_structure'] = table.index.isin(ids)
    table.to_csv(f'{W}/hits_annotated{suffix}.tsv', sep='\t')
    return table, explained, sim, ids


def scatter_pcoa(ax, t, color_col, colors, order, explained, legend_title, label_offset=(0.08, 0.12), legend_ncol=2):
    q = t.loc[QUERY]
    rest = t[t.has_structure].drop(QUERY)
    for k in order:
        s = rest[rest[color_col] == k]
        if s.empty:
            continue
        ax.scatter(s.PC1, s.PC2, s=4, lw=0, alpha=0.7, color=colors[k], label=f'{k} ({len(s)})', rasterized=True)
    ax.scatter(q.PC1, q.PC2, s=4, lw=0, color='black', zorder=10, label='Cinquedea (2090)')
    ax.annotate('cinquedea', (q.PC1, q.PC2), xytext=(q.PC1 + label_offset[0], q.PC2 + label_offset[1]), fontsize=7,
                ha='center', arrowprops={'arrowstyle': '-', 'lw': 0.6, 'color': 'black', 'shrinkB': 1})
    ax.set_xlabel(f'PCo1 ({explained[0]:.0%})')
    ax.set_ylabel(f'PCo2 ({explained[1]:.0%})')
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[::-1], labels[::-1], title=legend_title, title_fontsize=7, frameon=False, ncol=legend_ncol, fontsize=6.5,
                    loc='upper center', bbox_to_anchor=(0.5, -0.2), handletextpad=0.2, columnspacing=0.8)
    for h in leg.legendHandles:
        h.set_sizes([15])


def ranked_by_similarity(t):
    s = t.drop(QUERY).sort_values('tm_to_2090', ascending=False)
    s = s[s.has_structure].copy()
    s['similarity_rank'] = np.arange(1, len(s) + 1)
    return s


def panel_neighbours(ax, t, n=600):
    s = ranked_by_similarity(t)
    top = s[s.similarity_rank <= n]
    for k in DOMAIN_ORDER:
        d = top[top.domain == k]
        ax.scatter(d.similarity_rank, d.tm_to_2090, s=5, lw=0, color=DOMAIN_COLORS[k], label=k, rasterized=True)
    first_bact = s[s.domain == 'Bacteria'].iloc[0]
    ax.axvline(first_bact.similarity_rank, color='grey', lw=0.6, ls='--')
    ax.text(first_bact.similarity_rank + 8, 0.93, f'first bacterial hit\n(rank {first_bact.similarity_rank})',
            fontsize=7, va='top', color='grey')
    ax.set_xlabel('Rank of structural similarity to cinquedea')
    ax.set_ylabel('TM-score to cinquedea')
    return s


def figure_domain(t, explained, suffix='', label_offset=(0.08, 0.12)):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3), gridspec_kw={'width_ratios': [1, 1], 'wspace': 0.35})

    ax = axes[0]
    scatter_pcoa(ax, t, 'domain', DOMAIN_COLORS, DOMAIN_ORDER, explained, 'Domain', label_offset)
    ax.set_title('A  Structural similarity space', loc='left', fontweight='bold')

    ax = axes[1]
    s = panel_neighbours(ax, t)
    ax.set_title('B  Nearest structural neighbours', loc='left', fontweight='bold')

    fig.savefig(f'{FIG}/fig_foldseek_domain{suffix}.pdf', bbox_inches='tight')
    fig.savefig(f'{FIG}/fig_foldseek_domain{suffix}.png', bbox_inches='tight', dpi=300)
    return s


def figure_combined(t, explained, suffix='', label_offset=(0.08, 0.12)):
    """PCoA by domain, PCoA by Pfam domain and nearest structural neighbours, in one row."""
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), gridspec_kw={'wspace': 0.4})

    scatter_pcoa(axes[0], t, 'domain', DOMAIN_COLORS, DOMAIN_ORDER, explained, 'Domain', label_offset, legend_ncol=1)
    axes[0].set_title('A  Domain of life', loc='left', fontweight='bold')

    scatter_pcoa(axes[1], t, 'category', PFAM_COLORS, PFAM_ORDER[::-1], explained, 'Pfam domain', label_offset, legend_ncol=1)
    axes[1].set_title('B  Pfam domains', loc='left', fontweight='bold')

    panel_neighbours(axes[2], t)
    axes[2].set_title('C  Nearest structural neighbours', loc='left', fontweight='bold')

    fig.savefig(f'{FIG}/fig_foldseek_combined{suffix}.pdf', bbox_inches='tight')
    fig.savefig(f'{FIG}/fig_foldseek_combined{suffix}.png', bbox_inches='tight', dpi=300)


def figure_duf900(t, explained, suffix='', label_offset=(0.08, 0.12)):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3), gridspec_kw={'width_ratios': [1, 1.15], 'wspace': 0.35})

    ax = axes[0]
    scatter_pcoa(ax, t, 'category', PFAM_COLORS, PFAM_ORDER[::-1], explained, 'Pfam domain', label_offset)
    ax.set_title('A  Pfam domains in structural space', loc='left', fontweight='bold')

    ax = axes[1]
    s = t.drop(QUERY)
    s = s[s.has_structure]
    s = s[s.domain != 'Metagenome (unassigned)']
    s = s.assign(duf900=np.where(s.category == 'DUF900', 'DUF900', 'other'))
    domains = [d for d in ['Archaea', 'Bacteria', 'Eukaryota'] if d in set(s.domain)]
    groups = [(d, c) for d in domains for c in ['DUF900', 'other']]
    rng = np.random.default_rng(0)
    labels = []
    for i, (dom, cat) in enumerate(groups):
        v = s[(s.domain == dom) & (s.duf900 == cat)].tm_to_2090
        ax.scatter(i + rng.uniform(-0.3, 0.3, len(v)), v, s=3, lw=0, alpha=0.6,
                   color=DOMAIN_COLORS[dom], rasterized=True)
        ax.boxplot(v, positions=[i], widths=0.6, showfliers=False,
                   medianprops={'color': 'black'}, boxprops={'lw': 0.8}, whiskerprops={'lw': 0.8}, capprops={'lw': 0.8})
        labels.append(f'{cat}\n({len(v)})')
    for i, dom in enumerate(domains):
        ax.text(2 * i + 0.5, -0.2, dom, transform=ax.get_xaxis_transform(), ha='center', va='top',
                fontsize=7.5, color=DOMAIN_COLORS[dom], fontweight='bold')
    ax.set_xticks(range(len(groups)), labels, fontsize=6.5)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.6, len(groups) - 0.4)
    ax.set_ylabel('TM-score to cinquedea')
    ax.set_title('B  Similarity to cinquedea by DUF900 status', loc='left', fontweight='bold')

    fig.savefig(f'{FIG}/fig_foldseek_duf900{suffix}.pdf', bbox_inches='tight')
    fig.savefig(f'{FIG}/fig_foldseek_duf900{suffix}.png', bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    ava = load_all_vs_all()

    table, explained, sim, ids = build_table(ava)
    figure_domain(table, explained)
    figure_duf900(table, explained)
    figure_combined(table, explained)

    table, explained, sim, ids = build_table(ava, exclude=['Eukaryota'], suffix='_no_eukaryotes')
    figure_domain(table, explained, suffix='_no_eukaryotes', label_offset=(0.12, -0.08))
    figure_duf900(table, explained, suffix='_no_eukaryotes', label_offset=(0.12, -0.08))
    figure_combined(table, explained, suffix='_no_eukaryotes', label_offset=(0.12, -0.08))
