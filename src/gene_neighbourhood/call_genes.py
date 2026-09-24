"""Gene calling with pyrodigal (Prodigal single mode, translation table 11), as done by GTDB.

Gene identifiers follow the GTDB/Prodigal convention: <contig accession>_<gene number on contig>.

Run from the repo root:
    python src/gene_neighbourhood/call_genes.py
"""
import gzip
from multiprocessing import Pool
from pathlib import Path

import pandas as pd
import pyrodigal
from Bio import SeqIO

W = Path('data/gene_neighbourhood')


def call_genes(fna):
    genome = fna.name[:-len('.fna.gz')]
    with gzip.open(fna, 'rt') as f:
        contigs = [(r.id, bytes(r.seq)) for r in SeqIO.parse(f, 'fasta')]

    finder = pyrodigal.GeneFinder(meta=False)
    finder.train(*(s for _, s in contigs), translation_table=11)

    rows, faa = [], []
    for contig, seq in contigs:
        for i, gene in enumerate(finder.find_genes(seq), start=1):
            gene_id = f'{contig}_{i}'
            rows.append((genome, contig, len(seq), gene_id, i, gene.begin, gene.end,
                         '+' if gene.strand == 1 else '-', gene.partial_begin or gene.partial_end))
            faa.append(f'>{gene_id} {genome}\n{gene.translate(translation_table=11).rstrip("*")}\n')
    return rows, faa


if __name__ == '__main__':
    fnas = sorted((W / 'genomes').glob('*.fna.gz'))
    with Pool(8) as pool:
        results = pool.map(call_genes, fnas)

    genes = pd.DataFrame(
        [r for rows, _ in results for r in rows],
        columns=['genome', 'contig', 'contig_length', 'gene_id', 'gene_index', 'start', 'end', 'strand', 'partial'],
    )
    genes.to_csv(W / 'genes.tsv.gz', sep='\t', index=False)
    with gzip.open(W / 'proteins.faa.gz', 'wt') as f:
        for _, faa in results:
            f.writelines(faa)
    print(len(fnas), 'genomes', len(genes), 'genes')
