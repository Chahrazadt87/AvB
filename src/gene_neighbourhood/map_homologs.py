"""Locate each cinquedea homolog among the predicted genes of its own genome (by sequence).

GTDB gene numbers occasionally differ from ours (e.g. different assembly version), so homologs
are mapped by sequence search against the genes of their genome rather than by identifier.

Run from the repo root:
    python src/gene_neighbourhood/map_homologs.py
"""
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

W = Path('data/gene_neighbourhood')
MMSEQS = Path.home() / 'bin' / 'mmseqs'
QUERY_ID = 'pgaptmp_002090_1'
QUERY_GENOME = 's5a-1'


def homologs():
    s = pd.read_csv('data/supplementary/Table_S4_sequence_homologs.tsv', sep='\t')
    s.loc[s.protein_id == QUERY_ID, 'genome_accession'] = QUERY_GENOME
    return s[s.genome_accession.notna()].reset_index(drop=True)


if __name__ == '__main__':
    s = homologs()
    s['query'] = [f'q{i}' for i in range(len(s))]

    with tempfile.TemporaryDirectory() as tmp:
        faa = Path(tmp) / 'homologs.faa'
        faa.write_text(''.join(f'>{q}\n{seq.rstrip("*")}\n' for q, seq in zip(s['query'], s.sequence)))
        subprocess.run([
            MMSEQS, 'easy-search', faa, W / 'proteins.faa.gz', Path(tmp) / 'hits.tsv', Path(tmp) / 'tmp',
            '--threads', '8', '-s', '7.5', '--max-seqs', '2000',
            '--format-output', 'query,target,fident,alnlen,qlen,tlen,evalue,bits',
        ], check=True, capture_output=True)
        hits = pd.read_csv(Path(tmp) / 'hits.tsv', sep='\t', header=None,
                           names=['query', 'gene_id', 'fident', 'alnlen', 'qlen', 'tlen', 'evalue', 'bits'])

    genes = pd.read_csv(W / 'genes.tsv.gz', sep='\t', usecols=['genome', 'gene_id']).set_index('gene_id').genome
    hits['genome'] = hits.gene_id.map(genes)
    hits = hits.merge(s[['query', 'genome_accession', 'gtdb_gene_id']], on='query')
    # Identical copies can exist in a genome: on ties, prefer the gene reported by GTDB.
    hits['is_gtdb_gene'] = hits.gene_id == hits.gtdb_gene_id
    best = (hits[hits.genome == hits.genome_accession]
            .sort_values(['bits', 'is_gtdb_gene'], ascending=False)
            .drop_duplicates('query'))

    out = s[['query', 'protein_id', 'genome_accession', 'gtdb_gene_id', 'domain', 'gtdb_class', 'gtdb_species',
             'mmseqs_bitscore']].merge(best[['query', 'gene_id', 'fident', 'qlen', 'tlen']], on='query', how='left')
    out = out.drop(columns='query').rename(columns={'genome_accession': 'genome', 'mmseqs_bitscore': 'search_bitscore'})
    out.to_csv(W / 'homolog_genes.tsv', sep='\t', index=False)

    print('mapped', out.gene_id.notna().sum(), 'of', len(out), '| min identity', out.fident.min(),
          '| differs from GTDB id', (out.gene_id != out.gtdb_gene_id).sum() - 1)
