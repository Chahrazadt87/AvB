"""Gene neighbourhood of cinquedea in Halobacteria: gene maps only.

    (A) gene maps around the top 10 homologs, one per species
    (B) gene maps around the top 10 homologs, one per genus

Families shared by the displayed neighbourhoods have the same colour in A and B. Homologs whose contig
does not extend MIN_CONTIG_ROOM_KB on both sides are skipped, so that no gene map is cut by a contig edge.
Gene maps are ordered by the GTDB r214 archaeal tree (cladogram on the left), H. larsenii s5a-1 first.

Run from the repo root:
    python src/gene_neighbourhood/make_halobacteria_maps.py [window in kb, default 10]
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory

sys.path.insert(0, 'src/gene_neighbourhood')
from make_figures import (FIG, draw_gene_maps, family_legend_handles, family_names, gene_map_windows,
                          select_homologs, shared_family_colors)
from make_halobacteria_figure import CLADE
from neighbourhood import GTDB, QUERY_GENOME, load_family_pfam, load_genes, load_homologs, load_taxonomy

N_HOMOLOGS = 10
N_COLORS = 16
# Largest window, so that the 5, 10 and 20 kb figures show the same genomes
MIN_CONTIG_ROOM_KB = 20
# Colours given to other families, per window, as {gene of the family losing the colour: gene of the family
# gaining it} (genes of H. larsenii s5a-1). 10 kb (main figure): the conserved neighbours at -4
# (Glyco_transf_4) and -5 (Copper-bind) are coloured instead of MeaB + cobW (-12) and DnaJ (+8).
RECOLOR = {10: {'s5a1_contig_1_2087': 's5a1_contig_1_2079', 's5a1_contig_1_2067': 's5a1_contig_1_2080'}}


def select_complete(homologs, genes, by):
    """Top homologs in Halobacteria, one per species or genus, with a complete neighbourhood.

    E.g. excludes Halalkalicoccus sp022561895 (GCA_022561895.1), whose homolog is on a 7 kb contig.
    """
    g = genes.set_index('gene_id').join(genes.groupby('contig').agg(lo=('start', 'min'), hi=('end', 'max')), on='contig')
    room = np.minimum(g.start - g.lo, g.hi - g.end)
    complete = (room.reindex(homologs.gene_id) >= MIN_CONTIG_ROOM_KB * 1000).values
    return select_homologs(homologs[complete], by, N_HOMOLOGS, clade=CLADE)


def family_colors(windows, genes, window_kb):
    """Coloured families and their colours (see shared_family_colors), with the swaps of RECOLOR."""
    shared, colors = shared_family_colors(windows, N_COLORS)
    family = genes.set_index('gene_id').family
    for old, new in RECOLOR.get(window_kb, {}).items():
        old, new = family[old], family[new]
        shared[shared.index(old)] = new
        colors[new] = colors.pop(old)
    return shared, colors


def species_tree(top):
    """GTDB r214 tree of the selected genomes, and the genomes reordered by their position in the tree.

    Leaves are named by genome. H. larsenii s5a-1 is not in GTDB and is placed as a sister of its
    species representative (see Methods, GTDB-Tk). At each node, the clade with the best homolog comes
    first, so s5a-1 is on top and the other genomes follow by increasing phylogenetic distance from it.
    """
    from ete3 import Tree
    rep = load_taxonomy(top.genome).gtdb_genome_representative
    t = Tree(str(GTDB / 'ar53_r214.tree'), format=1, quoted_node_names=True)
    t.prune([t & r for r in set(rep)])
    for leaf in t.get_leaves():
        genomes = list(rep.index[rep == leaf.name])
        if len(genomes) == 1:
            leaf.name = genomes[0]
        else:
            leaf.name = ''
            for g in genomes:
                leaf.add_child(name=g)

    score = top.set_index('genome').search_bitscore
    for node in t.traverse('postorder'):
        node.add_feature('best', max(score[g] for g in node.get_leaf_names()))
        node.children.sort(key=lambda c: -c.best)
    assert t.get_leaf_names()[0] == QUERY_GENOME
    return t, top.set_index('genome').loc[t.get_leaf_names()].reset_index()


def draw_cladogram(ax, tree, tip_y, lw=0.6):
    """Cladogram with tips aligned on the right (x=0) at the given figure heights."""
    tr = blended_transform_factory(ax.transData, ax.figure.transFigure)
    pos = {}
    for node in tree.traverse('postorder'):
        if node.is_leaf():
            pos[node] = (0, tip_y[node.name])
            continue
        xs, ys = zip(*(pos[c] for c in node.children))
        x = min(xs) - 1
        pos[node] = (x, (ys[0] + ys[-1]) / 2)
        ax.plot([x, x], [ys[0], ys[-1]], color='black', lw=lw, transform=tr, solid_capstyle='butt')
        for cx, cy in zip(xs, ys):
            ax.plot([x, cx], [cy, cy], color='black', lw=lw, transform=tr, solid_capstyle='butt')
    x, y = pos[tree]
    ax.plot([x - 0.4, x], [y, y], color='black', lw=lw, transform=tr)
    ax.set_xlim(x - 0.4, 0)
    ax.axis('off')


if __name__ == '__main__':
    window_kb = int(sys.argv[1]) if len(sys.argv) > 1 else 10

    genes = load_genes()
    homologs = load_homologs()
    pfam = load_family_pfam()

    tree_species, top_species = species_tree(select_complete(homologs, genes, 'species'))
    tree_genus, top_genus = species_tree(select_complete(homologs, genes, 'genus'))
    by_species = gene_map_windows(genes, top_species, window_kb)
    by_genus = gene_map_windows(genes, top_genus, window_kb)
    shared, colors = family_colors(by_species + by_genus, genes, window_kb)

    fig = plt.figure(figsize=(7.2, 4.6))
    gs = fig.add_gridspec(1, 2, wspace=0.08)
    for col, (windows, tree, title) in enumerate([(by_species, tree_species, 'A  Top homologs, one per species'),
                                                  (by_genus, tree_genus, 'B  Top homologs, one per genus')]):
        panel = gs[0, col].subgridspec(1, 2, width_ratios=[0.09, 1], wspace=0.03)
        sub = panel[0, 1].subgridspec(N_HOMOLOGS, 1, hspace=0.9)
        axes = [fig.add_subplot(sub[i]) for i in range(N_HOMOLOGS)]
        draw_gene_maps(axes, windows, colors, window_kb, species_fontsize=6.5, genome_fontsize=5.5)
        # Tips at the height of the gene line of each map
        tip_y = {h.genome: ax.get_position().y0 + ax.get_position().height / 2 for ax, (h, _, _) in zip(axes, windows)}
        tree_ax = fig.add_subplot(panel[0, 0])
        draw_cladogram(tree_ax, tree, tip_y)
        map_box, tree_box = axes[0].get_position(), tree_ax.get_position()
        axes[0].set_title(title, loc='left', fontweight='bold', fontsize=8, pad=14,
                          x=(tree_box.x0 - map_box.x0) / map_box.width)

    handles = family_legend_handles(shared, colors, pfam, family_names(shared, pfam))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 0.06), ncol=3, frameon=False,
               fontsize=6.5, title=f'Gene families (±{window_kb} kb around the cinquedea homolog)', title_fontsize=7)

    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_maps_{window_kb}kb.pdf', bbox_inches='tight')
    fig.savefig(FIG / f'fig_neighbourhood_halobacteria_maps_{window_kb}kb.png', bbox_inches='tight', dpi=300)
