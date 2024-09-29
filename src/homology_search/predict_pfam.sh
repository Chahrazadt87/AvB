#!/bin/bash
set -e

PfamA_37=""  # Download from https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.hmm.gz
search_output="data/hydrolase_search/search_output.fasta"
hmmer_output="data/hydrolase_search/search_output.pfam.domtblout.txt"

hmmsearch \
    -o /dev/null \
    --domtblout $hmmer_output \
    --cut_ga \
    --cpu 4 \
    $PfamA_37
    $search_output
