#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -e error_mmseqs_db_search.txt
#PBS -o output_mmseqs_db_search.txt

set -e
start=$(date +%s%N)

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate mmseqs

# Input and outputs files are passed as arguments to qsub
if [ -z "${INPUT_FASTA}" ]; then
    echo "INPUT_FASTA variable is not set: use qsub -v to set it"
    exit 1
fi

if [ "${TARGET_DB}" = "UniProtKB" ]; then
    TARGET_DB=/rds/general/project/lms-warnecke-raw/live/UniProtKB/UniProtKB
elif [ "${TARGET_DB}" = "GTDB_214" ]; then
    TARGET_DB=/rds/general/project/lms-warnecke-raw/live/GTDB_214/GTDB_214
elif [ "${TARGET_DB}" = "db_proka" ]; then
    TARGET_DB=/rds/general/project/lms-warnecke-raw/live/db_prokaryotes/mmseqs_db/db_proka
elif [ -z "${TARGET_DB}" ]; then
    echo "TARGET_DB variable is not set: use qsub -v to set it"
    exit 1
fi

if [ -z "${OUTPUT_FOLDER}" ]; then
    echo "OUTPUT_FOLDER variable is not set: use qsub -v to set it"
    exit 1
fi

if [ -z "${SENSITIVITY}" ]; then
    # MMSeqs2 default
    SENSITIVITY="5.7"
fi

compute_elapsed_time() {
    local start=$1
    local end=$2
    local elapsed_ns=$((end - start))  # elapsed time in nanoseconds
    local elapsed_ms=$((elapsed_ns / 1000000))  # convert to milliseconds
    local ms=$((elapsed_ms % 1000))  # milliseconds
    local total_seconds=$((elapsed_ms / 1000))  # total seconds
    local hours=$((total_seconds / 3600))
    local minutes=$(( (total_seconds % 3600) / 60 ))
    local seconds=$((total_seconds % 60))

    local output=""

    if [ $hours -gt 0 ]; then
        output="${hours}h "
    fi

    if [ $minutes -gt 0 ] || [ $hours -gt 0 ]; then
        output="${output}${minutes}m "
    fi

    if [ $seconds -gt 0 ] || [ $minutes -gt 0 ] || [ $hours -gt 0 ]; then
        output="${output}${seconds}s "
    fi

    if [ $total_seconds -lt 1 ]; then
        output="${ms}ms"
    elif [ $ms -gt 0 ]; then
        output="${output}${ms}ms"
    fi

    echo $output
}

TMP=/rds/general/user/rs1521/ephemeral
TEMP_DIR=$(mktemp -d "${TMP}/tmp.XXXXXX")

echo "Searching for proteins in ${TARGET_DB} using mmseqs2"
echo "Input FASTA: ${INPUT_FASTA}"
echo "Output folder: ${OUTPUT_FOLDER}"
echo "Temp directory (will be deleted upon completion): ${TEMP_DIR}"

mmseqs createdb \
	${INPUT_FASTA} \
	${TEMP_DIR}/input_db

mmseqs search \
	${TEMP_DIR}/input_db \
	${TARGET_DB} \
	${TEMP_DIR}/results \
	${TMP} \
	--threads 32 --split-memory-limit "125G" -s "${SENSITIVITY}"

mmseqs convertalis \
	${TEMP_DIR}/input_db \
    ${TARGET_DB} \
    ${TEMP_DIR}/results \
	${OUTPUT_FOLDER}/results.tsv \
	--threads 32 --format-mode 4 --format-output "query,target,evalue,bits,tstart,tend,taxlineage"

# Export list of target ids
python - <<END
import pandas as pd

df = pd.read_csv("${OUTPUT_FOLDER}/results.tsv", sep="\t")
unique_values = sorted(df['target'].unique())

with open("${TEMP_DIR}/result_ids.txt", 'w') as f:
    for val in unique_values:
        f.write(str(val) + '\n')
END

mmseqs createsubdb \
	${TEMP_DIR}/result_ids.txt \
	${TARGET_DB} \
	${OUTPUT_FOLDER}/results_db \
	--id-mode 1

mmseqs convert2fasta \
	${OUTPUT_FOLDER}/results_db \
	${OUTPUT_FOLDER}/results.fasta

rm -r ${TEMP_DIR}

end=$(date +%s%N)
echo "Total elapsed time: $(compute_elapsed_time $start $end)"

exit 0
