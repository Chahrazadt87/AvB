"""Tables of the genes shown in the Halobacteria gene maps (fig_neighbourhood_halobacteria_maps_<N>kb).

One workbook per window, with the genes of panel A (top homologs, one per species) and panel B
(one per genus). Our Prodigal genes are matched to the NCBI annotation of the same assembly
(same contig and strand, same stop codon or >=50% overlap); NCBI features without a Prodigal match
(tRNAs, rRNAs, pseudogenes, ...) are listed too. Each Prodigal gene is annotated with Pfam (release 37).

Run from the repo root (after download_annotations.sh on the genomes of the gene maps):
    python src/gene_neighbourhood/export_gene_maps.py

Outputs (data/gene_neighbourhood/tables/):
    gene_maps_halobacteria_<N>kb.xlsx   for N = 5, 10, 20
"""
import gzip
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

import numpy as np
import pandas as pd
from Bio import SeqIO

sys.path.insert(0, 'src/gene_neighbourhood')
from make_figures import family_names, gene_map_windows
from make_halobacteria_maps import family_colors, select_complete, species_tree
from neighbourhood import QUERY_GENOME, W, load_family_pfam, load_genes, load_homologs

HMMSEARCH = '/Users/rs1521/miniconda3_x86/bin/hmmsearch'
PFAM_HMM = Path('/Users/rs1521/workspace/data/Pfam_37/Pfam-A.hmm')
ANNOTATIONS = W / 'annotations'
# PGAP annotation of H. larsenii s5a-1 (contig 1 only, named 'Haloferax' in the GFF)
QUERY_GFF = Path('data/gene_context/genomes/Haloferax_larsenii_S5a1/annot.gff')
PFAM_DOMTBL = W / 'gene_maps.pfam.domtblout.txt'
WINDOWS_KB = [5, 10, 20]
OUT = W / 'tables'

# Columns of the tables (internal name: column name). Genes without a locus tag have no NCBI match;
# rows without a Prodigal gene id are NCBI features not drawn in the figure (e.g. tRNAs).
COLUMNS = {
    'panel_row': 'panel_row',
    'species': 'species',
    'genome': 'genome',
    'annotation': 'annotation',
    'position': 'position_relative_to_homolog',
    'is_cinquedea_homolog': 'is_cinquedea_homolog',
    'contig': 'contig',
    'start': 'start',
    'end': 'end',
    'strand': 'strand',
    'gene_id': 'prodigal_gene_id',
    'length_aa': 'length_aa',
    'pfam_domains': 'pfam_domains',
    'ncbi_locus_tag': 'locus_tag',
    'ncbi_protein_id': 'protein_id',
    'ncbi_product': 'description',
}


# --- NCBI annotation ----------------------------------------------------------------------------

def read_gff(path, genome):
    """Features of an NCBI GFF3, one row per locus (gene coordinates, product feature attributes)."""
    rows = []
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else open(path)) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            x = line.rstrip('\n').split('\t')
            attr = dict(kv.split('=', 1) for kv in x[8].split(';') if '=' in kv)
            rows.append((x[0], x[2], int(x[3]), int(x[4]), x[6],
                         *(unquote(attr.get(k, '')) for k in ['locus_tag', 'old_locus_tag', 'protein_id',
                                                                 'product', 'Note', 'pseudo'])))
    g = pd.DataFrame(rows, columns=['contig', 'type', 'start', 'end', 'strand', 'locus_tag', 'old_locus_tag',
                                    'protein_id', 'product', 'note', 'pseudo'])
    g = g[g.locus_tag != '']
    genes = g[g.type.isin(['gene', 'pseudogene'])].drop_duplicates('locus_tag')
    products = (g[~g.type.isin(['gene', 'pseudogene', 'exon'])]
                .groupby('locus_tag').agg(ncbi_feature=('type', 'first'), protein_id=('protein_id', 'first'),
                                          product=('product', 'first'), note=('note', 'first'),
                                          pseudo=('pseudo', 'max')))
    genes = genes.drop(columns=['type', 'protein_id', 'product', 'note', 'pseudo']).merge(
        products, left_on='locus_tag', right_index=True, how='left')
    genes['pseudo'] = genes.pseudo.eq('true') | genes.locus_tag.isin(g.locus_tag[g.type == 'pseudogene'])
    genes['protein_id'] = genes.protein_id.str.replace('extdb:', '', regex=False)
    genes['genome'] = genome
    return genes.rename(columns={c: f'ncbi_{c}' for c in ['start', 'end', 'strand', 'locus_tag', 'old_locus_tag',
                                                           'protein_id', 'product', 'note', 'pseudo']})


def load_annotation(genome):
    if genome == QUERY_GENOME:
        a = read_gff(QUERY_GFF, genome)
        a['contig'] = a.contig.replace({'Haloferax': 's5a1_contig_1'})
        return a
    return read_gff(ANNOTATIONS / f'{genome}.gff.gz', genome)


def match_annotation(c, ncbi):
    """Match Prodigal genes of a window to NCBI features (stop codon, then overlap on the same strand)."""
    stop = lambda s, e, st: np.where(st == '+', e, s)
    c = c.assign(stop=stop(c.start, c.end, c.strand))
    n = ncbi.assign(stop=stop(ncbi.ncbi_start, ncbi.ncbi_end, ncbi.ncbi_strand))
    n = n[n.ncbi_feature.isin(['CDS']) | n.ncbi_pseudo].drop(columns='genome')

    m = c.merge(n, left_on=['contig', 'strand', 'stop'], right_on=['contig', 'ncbi_strand', 'stop'], how='left')
    m['ncbi_match'] = np.where(m.ncbi_locus_tag.isna(), 'none',
                               np.where((m.start == m.ncbi_start) & (m.end == m.ncbi_end), 'exact', 'same_stop'))
    m = m.drop_duplicates('gene_id')
    for i in m.index[m.ncbi_match == 'none']:
        x = m.loc[i]
        o = n[(n.contig == x.contig) & (n.ncbi_strand == x.strand)]
        ov = np.minimum(o.ncbi_end, x.end) - np.maximum(o.ncbi_start, x.start) + 1
        frac = ov / np.minimum(o.ncbi_end - o.ncbi_start + 1, x.end - x.start + 1)
        if len(o) and frac.max() >= 0.5:
            best = o.loc[frac.idxmax()]
            for col in n.columns.drop(['contig', 'stop']):
                m.at[i, col] = best[col]
            m.at[i, 'ncbi_match'] = 'overlap'
    return m.drop(columns='stop')


# --- Pfam ---------------------------------------------------------------------------------------

def run_pfam(gene_ids):
    if PFAM_DOMTBL.exists():
        return
    faa = W / 'gene_maps.faa'
    with gzip.open(W / 'proteins.faa.gz', 'rt') as f:
        SeqIO.write((r for r in SeqIO.parse(f, 'fasta') if r.id in gene_ids), faa, 'fasta')
    subprocess.run([HMMSEARCH, '-o', '/dev/null', '--domtblout', PFAM_DOMTBL, '--cut_ga', '--cpu', '8',
                    PFAM_HMM, faa], check=True)
    faa.unlink()


def gene_pfam():
    rows = []
    with open(PFAM_DOMTBL) as f:
        for line in f:
            if not line.startswith('#'):
                x = line.split()
                rows.append((x[0], x[3], int(x[19])))
    p = pd.DataFrame(rows, columns=['gene_id', 'name', 'start']).sort_values(['gene_id', 'start'])
    return p.groupby('gene_id').agg(pfam_domains=('name', lambda s: ' + '.join(s)))


# --- Tables -------------------------------------------------------------------------------------

def protein_lengths(gene_ids):
    with gzip.open(W / 'proteins.faa.gz', 'rt') as f:
        return {r.id: len(r.seq) for r in SeqIO.parse(f, 'fasta') if r.id in gene_ids}


def window_table(h, a, c, window_kb, ncbi, pfam, lengths, labels, row):
    m = match_annotation(c, ncbi)
    m['in_figure'] = True

    # NCBI features of the window without a Prodigal gene
    span = window_kb * 1000
    n = ncbi[(ncbi.contig == a.contig) & (ncbi.ncbi_end >= a.start - span) & (ncbi.ncbi_start <= a.end + span)]
    extra = n[~n.ncbi_locus_tag.isin(m.ncbi_locus_tag.dropna())].copy()
    extra = extra.assign(ncbi_match='ncbi_only', in_figure=False, start=extra.ncbi_start, end=extra.ncbi_end,
                         strand=extra.ncbi_strand)
    t = pd.concat([m, extra], ignore_index=True)

    flip = a.strand == '-'
    t = t.assign(key=np.where(flip, -t.end, t.start)).sort_values('key').reset_index(drop=True)
    t['is_cinquedea_homolog'] = t.gene_id == h.gene_id
    t['position'] = t.index - t.index[t.is_cinquedea_homolog][0]
    t['orientation'] = np.where(t.strand == a.strand, 'same', 'opposite')
    t['distance_to_homolog_bp'] = np.maximum(0, np.maximum(t.start - a.end, a.start - t.end))
    t['length_aa'] = t.gene_id.map(lengths)
    t['figure_label'] = np.where(t.is_cinquedea_homolog, 'Cinquedea homolog',
                                 t.family.map(labels).fillna('Other')).astype(object)
    t.loc[~t.in_figure, 'figure_label'] = ''
    t = t.join(pfam, on='gene_id')
    t['panel_row'] = row
    t['species'] = h.species
    t['genome'] = h.genome
    t['annotation'] = ('PGAP (this study)' if h.genome == QUERY_GENOME
                       else 'RefSeq' if h.genome.startswith('GCF_') else 'GenBank')
    t['partial'] = t.partial.astype('boolean')
    t['ncbi_pseudo'] = t.ncbi_pseudo.where(t.ncbi_locus_tag.notna()).astype('boolean')
    return t[list(COLUMNS)]


if __name__ == '__main__':
    genes = load_genes()
    homologs = load_homologs()
    family_pfam = load_family_pfam()

    # Same order as in the figure (GTDB tree)
    panels = {'A': species_tree(select_complete(homologs, genes, 'species'))[1],
              'B': species_tree(select_complete(homologs, genes, 'genus'))[1]}
    genomes = pd.concat(panels.values()).genome.unique()
    ncbi = {g: load_annotation(g) for g in genomes}

    widest = [gene_map_windows(genes, top, max(WINDOWS_KB)) for top in panels.values()]
    gene_ids = {gid for windows in widest for _, _, c in windows for gid in c.gene_id}
    run_pfam(gene_ids)
    pfam = gene_pfam()
    lengths = protein_lengths(gene_ids)

    OUT.mkdir(exist_ok=True)
    for window_kb in WINDOWS_KB:
        windows = {p: gene_map_windows(genes, top, window_kb) for p, top in panels.items()}
        # Same colours (legend entries) as in the figure
        shared, _ = family_colors(windows['A'] + windows['B'], genes, window_kb)
        labels = family_names(shared, family_pfam)
        tables = {p: pd.concat([window_table(h, a, c, window_kb, ncbi[h.genome], pfam, lengths, labels, i + 1)
                                for i, (h, a, c) in enumerate(w)], ignore_index=True)
                  for p, w in windows.items()}

        path = OUT / f'gene_maps_halobacteria_{window_kb}kb.xlsx'
        with pd.ExcelWriter(path) as w:
            tables['A'].rename(columns=COLUMNS).to_excel(w, sheet_name='A one per species', index=False)
            tables['B'].rename(columns=COLUMNS).to_excel(w, sheet_name='B one per genus', index=False)
            for sheet in w.sheets.values():
                sheet.freeze_panes = 'A2'
        for p, t in tables.items():
            print(path.name, p, len(t), 'rows |', t.ncbi_locus_tag.notna().sum(), 'with locus tag')
