#!/bin/bash
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -e error_hydrolase_2090_search.txt
#PBS -o output_hydrolase_2090_search.txt

set -e
start=$(date +%s)

cd $HOME
mkdir -p chahrazad

INPUT="${HOME}/chahrazad/2090.fa"
GTDB_OUTPUT="chahrazad/2090_search_gtdb"
UNIPROT_OUTPUT="chahrazad/2090_search_uniprot"
export SENSITIVITY="7.5"

mkdir -p $GTDB_OUTPUT
mkdir -p $UNIPROT_OUTPUT

echo "Search in GTDB"
export INPUT_FASTA=$INPUT
export TARGET_DB="GTDB_214"
export OUTPUT_FOLDER=$GTDB_OUTPUT
bash mmseqs_db_search.sh

echo "Search in UniProtKB"
export INPUT_FASTA=$INPUT
export TARGET_DB="UniProtKB"
export OUTPUT_FOLDER=$UNIPROT_OUTPUT
bash mmseqs_db_search.sh

echo "Map GTDB results to UniProtKB"
export INPUT_FASTA="${GTDB_OUTPUT}/results.fasta"
export TARGET_DB="UniProtKB"
export OUTPUT_TSV="${UNIPROT_OUTPUT}/results_uniprot_mapped_to_gtdb.tsv"
bash map_sequences_to_db.sh

echo "Map GTDB results to our 'db_proka' database (Strock et al., 2024 - DOI: 10.1101/2024.09.18.613068)"
export INPUT_FASTA="${GTDB_OUTPUT}/results.fasta"
export TARGET_DB="db_proka"
export OUTPUT_TSV="${GTDB_OUTPUT}/results_gtdb_mapped_to_db_proka.tsv"
bash map_sequences_to_db.sh

end=$(date +%s)
echo "Total elapsed Time: $(($end-$start)) seconds"
