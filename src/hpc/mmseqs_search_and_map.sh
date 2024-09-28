#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -e error_mmseqs_search_and_map_gtdb_expanded.txt
#PBS -o output_mmseqs_search_and_map_gtdb_expanded.txt

set -e
start=$(date +%s)

cd /rds/general/user/rs1521/home/

export INPUT_FASTA="/rds/general/user/rs1521/home/chahrazad/2090_homologs_first_pass.fasta"
export TARGET_DB="GTDB_214"
export OUTPUT_FOLDER="/rds/general/user/rs1521/home/chahrazad/2090_search_gtdb_expanded"
export SENSITIVITY="7.5"

bash mmseqs_db_search.sh

export INPUT_FASTA="${OUTPUT_FOLDER}/results.fasta"
export TARGET_DB="db_proka"
export OUTPUT_TSV="${OUTPUT_FOLDER}/results_mapped_to_db_proka.tsv"

bash map_sequences_to_db.sh

end=$(date +%s)
echo "Total elapsed Time: $(($end-$start)) seconds"
