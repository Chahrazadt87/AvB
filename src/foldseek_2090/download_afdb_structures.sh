#!/bin/bash
# Download AFDB models for every Foldseek hit of hydrolase 2090.
# The October 2024 search used AFDB v4; the EBI now only serves v6 models.
set -e

hits="data/hydrolase_2090_AF3_foldseek_results.tsv"
outdir="data/foldseek_2090/structures"
mkdir -p $outdir

tail -n +2 $hits | cut -f2 | sed 's/-model_v4//' | sort -u | \
    xargs -P 16 -I{} sh -c \
    '[ -s '$outdir'/{}.pdb ] || curl -sf -o '$outdir'/{}.pdb https://alphafold.ebi.ac.uk/files/{}-model_v6.pdb || echo "missing {}"'
