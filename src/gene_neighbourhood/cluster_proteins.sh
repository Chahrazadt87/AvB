#!/bin/bash
# Cluster all predicted proteins of the 280 genomes into gene families.
set -e

W=data/gene_neighbourhood
TMP=$(mktemp -d)

~/bin/mmseqs easy-cluster \
    $W/proteins.faa.gz \
    $W/families \
    $TMP \
    --min-seq-id 0.3 \
    -c 0.8 \
    --cov-mode 0 \
    --threads 8

gzip -f $W/families_cluster.tsv
rm -f $W/families_all_seqs.fasta
rm -rf $TMP
