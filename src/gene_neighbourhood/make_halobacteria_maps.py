"""Gene neighbourhood of cinquedea in Halobacteria: gene maps only.

    (A) gene maps around the top 10 homologs, one per species
    (B) gene maps around the top 10 homologs, one per genus

Families shared by the displayed neighbourhoods have the same colour in A and B.

Run from the repo root:
    python src/gene_neighbourhood/make_halobacteria_maps.py [window in kb, default 10]
"""
import sys

import matplotlib.pyplot as plt

sys.path.insert(0, 'src/gene_neighbourhood')
from make_figures import (FIG, draw_gene_maps, family_legend_handles, family_names, gene_map_windows,
                          select_homologs, shared_family_colors)
from make_halobacteria_figure import CLADE
from neighbourhood import load_family_pfam, load_genes, load_homologs

N_HOMOLOGS = 10
N_COLORS = 16


if __name__ == '__main__':
    window_kb = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    genes = load_genes()
    homologs = load_homologs()
    pfam = load_family_pfam()

    by_species = gene_map_windows(genes, select_homologs(homologs, 'species', N_HOMOLOGS, clade=CLADE), window_kb)
    by_genus = gene_map_windows(genes, select_homologs(homologs, 'genus', N_HOMOLOGS, clade=CLADE), window_kb)
    shared, colors = shared_family_colors(by_species + by_genus, N_COLORS)

    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(1, 2, wspace=0.12)
    for col, (windows, title) in enumerate([(by_species, 'A  Top homologs, one per species'),
                                            (by_genus, 'B  Top homologs, one per genus')]):
        sub = gs[0, col].subgridspec(N_HOMOLOGS, 1, hspace=0.9)
        axes = [fig.add_subplot(sub[i]) for i in range(N_HOMOLOGS)]
        draw_gene_maps(axes, windows, colors, window_kb, species_fontsize=6.5, genome_fontsize=5.5)
        axes[0].set_title(title, loc='left', fontweight='bold', fontsize=8, pad=14)

    handles = family_legend_handles(shared, colors, pfam, family_names(shared, pfam))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.06), ncol=3, frameon=False,
               fontsize=6.5, title=f'Gene families (±{window_kb} kb around the cinquedea homolog)', title_fontsize=7)

    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_maps_{window_kb}kb.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_maps_{window_kb}kb.png', bbox_inches='tight', dpi=300)
