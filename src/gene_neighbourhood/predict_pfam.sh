#!/bin/bash
# Pfam (release 37) annotation of the representatives of families found in cinquedea neighbourhoods.
set -e

PfamA_37=/Users/rs1521/workspace/data/Pfam_37/Pfam-A.hmm
W=data/gene_neighbourhood

/opt/homebrew/Caskroom/miniforge/base/envs/amp/bin/python - <<'PY'
import pandas as pd
from Bio import SeqIO
W = 'data/gene_neighbourhood'
fams = set(pd.read_csv(f'{W}/window_20kb/neighbourhood_genes.tsv.gz', sep='\t').family.dropna())
fams |= set(pd.read_csv(f'{W}/window_20kb/control_families.tsv', sep='\t').family)
SeqIO.write((r for r in SeqIO.parse(f'{W}/families_rep_seq.fasta', 'fasta') if r.id in fams),
            f'{W}/neighbourhood_families.faa', 'fasta')
PY

/Users/rs1521/miniconda3_x86/bin/hmmsearch \
    -o /dev/null \
    --domtblout $W/neighbourhood_families.pfam.domtblout.txt \
    --cut_ga \
    --cpu 8 \
    $PfamA_37 \
    $W/neighbourhood_families.faa
