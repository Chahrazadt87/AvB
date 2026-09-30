#!/bin/bash
# Download NCBI genome annotations (GFF3) from NCBI Datasets, for the genomes shown in the gene maps.
set -e

genomes=$1   # one assembly accession per line
export outdir=data/gene_neighbourhood/annotations
mkdir -p $outdir

download() {
    acc=$1
    [ -s $outdir/$acc.gff.gz ] && return 0
    tmp=$(mktemp -d)
    if curl -sf -o $tmp/g.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_GFF" \
        && unzip -q -o $tmp/g.zip -d $tmp && [ -s $tmp/ncbi_dataset/data/$acc/genomic.gff ]; then
        gzip -c $tmp/ncbi_dataset/data/$acc/genomic.gff > $outdir/$acc.gff.gz
    else
        echo "failed $acc"
    fi
    rm -rf $tmp
}
export -f download

cat $genomes | xargs -P ${PARALLEL:-6} -n 1 bash -c 'download "$0"'
