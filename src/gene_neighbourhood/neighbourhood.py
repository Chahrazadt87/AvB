"""Gene neighbourhood analysis of cinquedea (hydrolase 2090) homologs.

Neighbourhood: all genes on the same contig within WINDOW bp of the anchor gene.
Gene families: MMseqs2 clusters (30% identity, 80% coverage) of all proteins of the genomes.
"""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

W = Path('data/gene_neighbourhood')
GTDB = Path('/Users/rs1521/workspace/data/gtdb_r214.1')
WINDOW = 20_000
QUERY_GENOME = 's5a-1'
SEED = 0


# --- Data ---------------------------------------------------------------------------------------

def load_genes():
    genes = pd.read_csv(W / 'genes.tsv.gz', sep='\t')
    fam = pd.read_csv(W / 'families_cluster.tsv.gz', sep='\t', header=None, names=['family', 'gene_id'])
    return genes.merge(fam, on='gene_id', how='left')


def load_homologs():
    return pd.read_csv(W / 'homolog_genes.tsv', sep='\t')


def load_taxonomy(genomes):
    """GTDB r214 taxonomy and species representative of each genome."""
    m = pd.concat([
        pd.read_csv(GTDB / f, sep='\t', low_memory=False,
                    usecols=['accession', 'gtdb_taxonomy', 'gtdb_genome_representative'])
        for f in ['ar53_metadata_r214.tsv', 'bac120_metadata_r214.tsv']
    ])
    m['genome'] = m.accession.str.replace(r'^(RS|GB)_', '', regex=True)
    ranks = m.gtdb_taxonomy.str.split(';', expand=True)
    for i, r in enumerate(['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']):
        m[r] = ranks[i].str[3:]
    m = m.set_index('genome').reindex(genomes)
    # H. larsenii s5a-1 is not in GTDB; it is placed in H. larsenii (see Methods, GTDB-Tk).
    if QUERY_GENOME in m.index:
        ref = pd.read_csv(GTDB / 'ar53_metadata_r214.tsv', sep='\t', low_memory=False,
                          usecols=['accession', 'gtdb_taxonomy', 'gtdb_genome_representative'])
        ref = ref[ref.gtdb_taxonomy.str.endswith('s__Haloferax larsenii')].iloc[0]
        ranks = ref.gtdb_taxonomy.split(';')
        m.loc[QUERY_GENOME, ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']] = [r[3:] for r in ranks]
        m.loc[QUERY_GENOME, 'gtdb_genome_representative'] = ref.gtdb_genome_representative
    m['clade'] = np.where(m['class'].isin(['Halobacteria', 'Nitrososphaeria']), m['class'],
                          np.where(m['domain'] == 'Archaea', 'Other archaea', m['domain']))
    return m[['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'clade', 'gtdb_genome_representative']]


# --- Neighbourhoods -----------------------------------------------------------------------------

def neighbourhoods(genes, anchor_ids, window=WINDOW):
    """Genes within `window` bp of each anchor (anchor excluded), with contig-edge truncation flags."""
    g = genes.set_index('gene_id')
    by_contig = {c: d for c, d in genes.groupby('contig')}
    rows, meta = [], []
    for a in anchor_ids:
        x = g.loc[a]
        lo, hi = x.start - window, x.end + window
        d = by_contig[x.contig]
        d = d[(d.end >= lo) & (d.start <= hi) & (d.gene_id != a)]
        rows.append(d.assign(anchor=a, offset=np.where(x.strand == '+', d.start - x.start, x.end - d.end),
                             same_strand=d.strand == x.strand))
        meta.append((a, x.genome, x.family, lo < 1 or hi > x.contig_length, len(d)))
    nb = pd.concat(rows, ignore_index=True)
    meta = pd.DataFrame(meta, columns=['anchor', 'genome', 'anchor_family', 'truncated', 'n_genes'])
    return nb, meta


def family_sets(nb):
    return nb.groupby('anchor').family.agg(lambda s: frozenset(s.dropna()))


# --- Anchors ------------------------------------------------------------------------------------

def cinquedea_anchors(homologs):
    """One anchor per genome: the homolog with the highest search score (cinquedea itself in s5a-1)."""
    h = homologs.dropna(subset=['gene_id']).copy()
    h['search_bitscore'] = h.search_bitscore.fillna(np.inf)
    return h.sort_values('search_bitscore', ascending=False).drop_duplicates('genome')


def random_anchors(genes, homolog_gene_ids, n_per_genome=20, seed=SEED):
    """Random genes (not homologs) used as background windows for the enrichment test."""
    pool = genes[~genes.gene_id.isin(homolog_gene_ids)]
    return pool.groupby('genome', group_keys=False).apply(
        lambda d: d.sample(min(n_per_genome, len(d)), random_state=seed))


def control_families(genes, taxonomy, cinquedea_families, n=10, min_prevalence=0.5, seed=SEED):
    """Random single-copy gene families spanning the same clades as cinquedea.

    Families must be present in a single copy in at least `min_prevalence` of the genomes of both
    Halobacteria and Nitrososphaeria.
    """
    g = genes.dropna(subset=['family'])
    copies = g.groupby(['family', 'genome']).size().rename('copies').reset_index()
    copies['clade'] = copies.genome.map(taxonomy.clade)
    n_genomes = taxonomy.clade.value_counts()
    single = copies[copies.copies == 1].groupby(['family', 'clade']).genome.nunique().unstack(fill_value=0)
    ok = pd.Series(True, index=single.index)
    for c in ['Halobacteria', 'Nitrososphaeria']:
        ok &= single.get(c, 0) >= min_prevalence * n_genomes[c]
    candidates = sorted(set(single.index[ok]) - set(cinquedea_families))
    rng = np.random.default_rng(seed)
    return list(rng.choice(candidates, size=min(n, len(candidates)), replace=False)), len(candidates)


def family_anchors(genes, family):
    """Members of a family in genomes where it is single-copy."""
    d = genes[genes.family == family]
    return d[~d.genome.duplicated(keep=False)]


# --- Phylogenetic distances ---------------------------------------------------------------------

def patristic_distances(taxonomy):
    """Patristic distances between species representatives in the GTDB r214 trees (per domain)."""
    from ete3 import Tree
    out = {}
    for domain, tree_file in [('Archaea', 'ar53_r214.tree'), ('Bacteria', 'bac120_r214.tree')]:
        reps = sorted(set(taxonomy[taxonomy.domain == domain].gtdb_genome_representative.dropna()))
        t = Tree(str(GTDB / tree_file), format=1, quoted_node_names=True)
        leaves = {l.name: l for l in t.iter_leaves() if l.name in reps}
        t.prune(list(leaves.values()), preserve_branch_length=True)
        leaves = {l.name: l for l in t.iter_leaves()}
        for a, b in combinations(sorted(leaves), 2):
            d = t.get_distance(leaves[a], leaves[b])
            out[(a, b)] = out[(b, a)] = d
    return out


def pairwise_conservation(sets, anchor_genome, taxonomy, distances, truncated):
    """Jaccard similarity of neighbourhood family sets for all pairs of genomes (same domain)."""
    rep = taxonomy.gtdb_genome_representative
    rows = []
    anchors = list(sets.index)
    for a, b in combinations(anchors, 2):
        ga, gb = anchor_genome[a], anchor_genome[b]
        ra, rb = rep.get(ga), rep.get(gb)
        if taxonomy.domain.get(ga) != taxonomy.domain.get(gb) or pd.isna(ra) or pd.isna(rb):
            continue
        d = 0.0 if ra == rb else distances.get((ra, rb))
        if d is None:
            continue
        sa, sb = sets[a], sets[b]
        union = len(sa | sb)
        rows.append((ga, gb, taxonomy.clade.get(ga), taxonomy.clade.get(gb), d,
                     len(sa & sb) / union if union else np.nan, truncated[a] or truncated[b]))
    return pd.DataFrame(rows, columns=['genome_a', 'genome_b', 'clade_a', 'clade_b', 'distance', 'jaccard', 'truncated'])


def load_family_pfam():
    """Pfam (release 37) domain architecture of each family representative, e.g. 'DnaJ+DnaJ_C'."""
    rows = []
    with open(W / 'neighbourhood_families.pfam.domtblout.txt') as f:
        for line in f:
            if not line.startswith('#'):
                x = line.split()
                rows.append((x[0], x[3], int(x[19])))
    p = pd.DataFrame(rows, columns=['family', 'pfam', 'start']).sort_values(['family', 'start'])
    return p.groupby('family').pfam.agg(lambda s: '+'.join(dict.fromkeys(s)))
