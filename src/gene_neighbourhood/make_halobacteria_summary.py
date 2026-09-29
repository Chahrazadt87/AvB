"""Summary figure of the gene neighbourhood of cinquedea in Halobacteria.

    (A) gene maps around the top 10 homologs, one per species
    (B) gene maps around the top 10 homologs, one per genus
    (C) cinquedea's ten most shared neighbour families (fraction of genomes), against the neighbours of the same
        rank of control gene families and random genes

Families shared by the displayed neighbourhoods have the same colour in A and B.

Run from the repo root:
    python src/gene_neighbourhood/make_halobacteria_summary.py [window in kb, default 10]
"""
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, 'src/gene_neighbourhood')
from make_figures import (FIG, draw_gene_maps, family_legend_handles, family_names, gene_map_windows,
                          select_homologs, shared_family_colors)
from make_halobacteria_figure import CLADE, compute_curves, plot_neighbour_bars
from neighbourhood import load_family_pfam, load_genes, load_homologs, load_taxonomy

N_HOMOLOGS = 10
N_COLORS = 16


if __name__ == '__main__':
    window_kb = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    genes = load_genes()
    homologs = load_homologs()
    taxonomy = load_taxonomy(sorted(genes.genome.unique()))
    pfam = load_family_pfam()

    by_species = gene_map_windows(genes, select_homologs(homologs, 'species', N_HOMOLOGS, clade=CLADE), window_kb)
    by_genus = gene_map_windows(genes, select_homologs(homologs, 'genus', N_HOMOLOGS, clade=CLADE), window_kb)
    shared, colors = shared_family_colors(by_species + by_genus, N_COLORS)
    curves = compute_curves(genes, homologs, taxonomy, window_kb * 1000)

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.55, 1], hspace=0.25, wspace=0.12, width_ratios=[1, 1])
    for col, (windows, title) in enumerate([(by_species, 'A  Top homologs, one per species'),
                                            (by_genus, 'B  Top homologs, one per genus')]):
        sub = gs[0, col].subgridspec(N_HOMOLOGS, 1, hspace=0.9)
        axes = [fig.add_subplot(sub[i]) for i in range(N_HOMOLOGS)]
        draw_gene_maps(axes, windows, colors, window_kb, species_fontsize=6.5, genome_fontsize=5.5)
        axes[0].set_title(title, loc='left', fontweight='bold', fontsize=8, pad=14)

    names = family_names(shared + list(curves['cinquedea_neighbours'].index[:10]), pfam)
    ax = fig.add_subplot(gs[1, 0])
    plot_neighbour_bars(ax, curves, window_kb, names, colors)
    ax.set_title(f'C  Most frequent neighbours of cinquedea (n={curves["n_genomes"]} genomes)', loc='left',
                 fontweight='bold', fontsize=8, x=-0.45)

    legend_ax = fig.add_subplot(gs[1, 1])
    legend_ax.axis('off')
    leg = legend_ax.legend(handles=family_legend_handles(shared, colors, pfam, names), loc='center left', frameon=False,
                           fontsize=6.5, ncol=1, title=f'Gene families in A and B (±{window_kb} kb)', title_fontsize=7)
    leg._legend_box.align = 'left'

    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_summary_{window_kb}kb.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_summary_{window_kb}kb.png', bbox_inches='tight', dpi=300)
