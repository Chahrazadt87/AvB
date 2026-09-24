#!/bin/bash
# Pfam annotation of all Foldseek hits (AFDB v6 sequences + archived sequences of retired UniProt entries).
set -e

PfamA_38=/Users/rs1521/workspace/data/Pfam_38/Pfam-A.hmm
W=data/foldseek_2090

cat $W/hits.fasta $W/retired.fasta > $W/all_hits.fasta

/Users/rs1521/miniconda3_x86/bin/hmmsearch \
    -o /dev/null \
    --domtblout $W/all_hits.pfam.domtblout.txt \
    --cut_ga \
    --cpu 4 \
    $PfamA_38 \
    $W/all_hits.fasta
