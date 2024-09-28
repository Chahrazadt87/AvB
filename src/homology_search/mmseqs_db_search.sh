#!/bin/bash
set -e

start=$(date +%s%N)

cd $HOME

. load_conda.sh
conda activate mmseqs

# Input and outputs files are passed as arguments to qsub
if [ -z "${INPUT_FASTA}" ]; then
    echo "INPUT_FASTA variable is not set: use qsub -v to set it"
    exit 1
fi

SHARED_DATA_DIR="/rds/general/project/lms-warnecke-raw/live"

if [ "${TARGET_DB}" = "UniProtKB" ]; then
    TARGET_DB="${SHARED_DATA_DIR}/UniProtKB/UniProtKB"
elif [ "${TARGET_DB}" = "GTDB_214" ]; then
    TARGET_DB="${SHARED_DATA_DIR}/GTDB_214/GTDB_214"
elif [ "${TARGET_DB}" = "db_proka" ]; then
    TARGET_DB="${SHARED_DATA_DIR}/db_prokaryotes/mmseqs_db/db_proka"
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

TEMP_DIR=$(mktemp -d "${HOME}/ephemeral/tmp.XXXXXX")

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
	${TEMP_DIR} \
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

end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"

exit 0
