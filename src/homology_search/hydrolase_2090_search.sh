#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -e error_hydrolase_2090_search.txt
#PBS -o output_hydrolase_2090_search.txt

set -e
start=$(date +%s)

cd $HOME
mkdir -p chahrazad
export SENSITIVITY="7.5"

# First pass
INPUT="${HOME}/chahrazad/2090.fa"
GTDB_OUTPUT="chahrazad/2090_search_gtdb"
UNIPROT_OUTPUT="chahrazad/2090_search_uniprot"

# Second pass
# INPUT="${HOME}/chahrazad/2090_homologs_first_pass.fasta"
# GTDB_OUTPUT = "chahrazad/2090_search_gtdb_expanded"
# UNIPROT_OUTPUT = "chahrazad/2090_search_uniprot_expanded"

mkdir -p $GTDB_OUTPUT
mkdir -p $UNIPROT_OUTPUT

export INPUT_FASTA=$INPUT
export TARGET_DB="GTDB_214"
export OUTPUT_FOLDER=$GTDB_OUTPUT
bash mmseqs_db_search.sh

export INPUT_FASTA="${GTDB_OUTPUT}/results.fasta"
export TARGET_DB="db_proka"
export OUTPUT_TSV="${GTDB_OUTPUT}/results_gtdb_mapped_to_db_proka.tsv"
bash map_sequences_to_db.sh

export INPUT_FASTA=$INPUT
export TARGET_DB="UniProtKB"
export OUTPUT_FOLDER=$UNIPROT_OUTPUT
bash mmseqs_db_search.sh

export INPUT_FASTA="${UNIPROT_OUTPUT}/results.fasta"
export TARGET_DB="db_proka"
export OUTPUT_TSV="${UNIPROT_OUTPUT}/results_uniprot_mapped_to_db_proka.tsv"
bash map_sequences_to_db.sh

export INPUT_FASTA="${GTDB_OUTPUT}/results.fasta"
export TARGET_DB="UniProtKB"
export OUTPUT_TSV="${UNIPROT_OUTPUT}/results_uniprot_mapped_to_gtdb.tsv"
bash map_sequences_to_db.sh

end=$(date +%s)
echo "Total elapsed Time: $(($end-$start)) seconds"
