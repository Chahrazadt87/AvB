#!/bin/bash
# All-vs-all Foldseek comparison of the hydrolase 2090 hits (+ the query itself).
set -e

FS=~/bin/foldseek/bin/foldseek  # version 463739e
W=data/foldseek_2090
TMP=$(mktemp -d)

$FS createdb $W/structures $W/db/hits --threads 8
$FS convert2fasta $W/db/hits $W/hits.fasta

$FS search $W/db/hits $W/db/hits $W/db/allvsall $TMP \
    --exhaustive-search 1 \
    -e inf \
    -a \
    --threads 8
$FS convertalis $W/db/hits $W/db/hits $W/db/allvsall $W/all_vs_all.tsv \
    --format-output query,target,fident,alnlen,qlen,tlen,evalue,bits,qtmscore,ttmscore,alntmscore

rm -rf $TMP
