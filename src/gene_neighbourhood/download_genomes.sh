#!/bin/bash
# Download genome sequences (FASTA) from NCBI Datasets for the genomes encoding cinquedea homologs.
set -e

genomes=$1   # one assembly accession per line
export outdir=data/gene_neighbourhood/genomes
mkdir -p $outdir

download() {
    acc=$1
    [ -s $outdir/$acc.fna.gz ] && return 0
    tmp=$(mktemp -d)
    if curl -sf -o $tmp/g.zip "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/$acc/download?include_annotation_type=GENOME_FASTA" \
        && unzip -q -o $tmp/g.zip -d $tmp; then
        cat $tmp/ncbi_dataset/data/$acc/*.fna | gzip > $outdir/$acc.fna.gz
    else
        echo "failed $acc"
    fi
    rm -rf $tmp
}
export -f download

cat $genomes | xargs -P ${PARALLEL:-6} -n 1 bash -c 'download "$0"'
