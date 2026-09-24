"""Gene neighbourhood analysis of cinquedea homologs: co-occurrence, enrichment and conservation.

Run from the repo root (after call_genes.py, map_homologs.py and cluster_proteins.sh):
    python src/gene_neighbourhood/run_analysis.py [window in kb, default 20]

Outputs (data/gene_neighbourhood/window_<N>kb/):
    neighbourhood_genes.tsv.gz   genes within ±N kb of each cinquedea homolog
    anchors.tsv                  one cinquedea homolog per genome, with taxonomy and window info
    family_enrichment.tsv        per clade: prevalence near cinquedea vs random windows (same genomes)
    conservation_pairs.tsv.gz    pairwise neighbourhood similarity vs phylogenetic distance
                                 (cinquedea and control families)
    control_families.tsv         control families
"""
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, 'src/gene_neighbourhood')
from neighbourhood import (W, WINDOW, cinquedea_anchors, control_families, family_anchors, family_sets,
                           load_genes, load_homologs, load_taxonomy, neighbourhoods,
                           pairwise_conservation, patristic_distances, random_anchors)

N_RANDOM_PER_GENOME = 20
N_CONTROL_FAMILIES = 10


def enrichment(nb, meta, bg_nb, bg_meta, taxonomy):
    rows = []
    for clade, m in meta.groupby(meta.genome.map(taxonomy.clade)):
        n = len(m)
        obs = nb[nb.anchor.isin(m.anchor)].drop_duplicates(['anchor', 'family']).family.value_counts()
        bgm = bg_meta[bg_meta.genome.isin(m.genome)]
        bg = bg_nb[bg_nb.anchor.isin(bgm.anchor)].drop_duplicates(['anchor', 'family']).family.value_counts()
        for fam, k in obs.items():
            p0 = (bg.get(fam, 0) + 0.5) / (len(bgm) + 1)
            rows.append((clade, fam, k, n, k / n, bg.get(fam, 0), len(bgm), p0,
                         binomtest(int(k), n, min(p0, 1), alternative='greater').pvalue))
    e = pd.DataFrame(rows, columns=['clade', 'family', 'n_windows_with_family', 'n_windows', 'prevalence',
                                    'n_random_windows_with_family', 'n_random_windows', 'random_prevalence', 'p'])
    e['fold_enrichment'] = e.prevalence / e.random_prevalence
    e['q'] = e.groupby('clade').p.transform(lambda p: multipletests(p, method='fdr_bh')[1])
    return e.sort_values(['clade', 'prevalence'], ascending=[True, False])


if __name__ == '__main__':
    window = int(sys.argv[1]) * 1000 if len(sys.argv) > 1 else WINDOW
    out = W / f'window_{window // 1000}kb'
    out.mkdir(exist_ok=True)

    genes = load_genes()
    homologs = load_homologs()
    taxonomy = load_taxonomy(sorted(genes.genome.unique()))
    cinquedea_families = set(genes.set_index('gene_id').loc[homologs.gene_id, 'family'])

    # Cinquedea neighbourhoods
    anchors = cinquedea_anchors(homologs)
    nb, meta = neighbourhoods(genes, anchors.gene_id, window)
    nb['is_cinquedea_family'] = nb.family.isin(cinquedea_families)
    nb.to_csv(out / 'neighbourhood_genes.tsv.gz', sep='\t', index=False)
    meta = meta.merge(anchors[['gene_id', 'protein_id', 'search_bitscore']], left_on='anchor', right_on='gene_id').drop(columns='gene_id')
    meta = meta.join(taxonomy, on='genome')
    meta.to_csv(out / 'anchors.tsv', sep='\t', index=False)

    # Enrichment against random windows in the same genomes
    bg_anchors = random_anchors(genes, set(homologs.gene_id), n_per_genome=N_RANDOM_PER_GENOME)
    bg_nb, bg_meta = neighbourhoods(genes, bg_anchors.gene_id, window)
    enrichment(nb, meta, bg_nb, bg_meta, taxonomy).to_csv(out / 'family_enrichment.tsv', sep='\t', index=False)

    # Conservation vs phylogenetic distance: cinquedea and control families
    distances = patristic_distances(taxonomy)
    controls, n_candidates = control_families(genes, taxonomy, cinquedea_families, n=N_CONTROL_FAMILIES)
    pairs = []
    for name, anchor_ids in [('cinquedea', meta.anchor)] + [(f, family_anchors(genes, f).gene_id) for f in controls]:
        f_nb, f_meta = (nb, meta) if name == 'cinquedea' else neighbourhoods(genes, anchor_ids, window)
        sets = family_sets(f_nb).reindex(f_meta.anchor).apply(lambda s: s if isinstance(s, frozenset) else frozenset())
        p = pairwise_conservation(sets, f_meta.set_index('anchor').genome, taxonomy, distances,
                                  f_meta.set_index('anchor').truncated)
        pairs.append(p.assign(anchor_family=name))
    pd.concat(pairs).to_csv(out / 'conservation_pairs.tsv.gz', sep='\t', index=False)

    pd.DataFrame({'family': controls}).to_csv(out / 'control_families.tsv', sep='\t', index=False)
    print(f'±{window // 1000} kb: {len(meta)} anchors; {len(bg_meta)} random windows; {len(controls)} control families '
          f'(from {n_candidates} candidates)')
