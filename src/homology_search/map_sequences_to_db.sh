#!/bin/bash
set -e

cd $HOME

. load_conda.sh
conda activate mmseqs

start=$(date +%s)

##
# Parameters are passed as arguments to qsub
##
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

if [ -z "${OUTPUT_TSV}" ]; then
    echo "OUTPUT variable is not set: use qsub -v to set it"
    exit 1
fi

TMP=$(mktemp -d "${HOME}/ephemeral/tmp.XXXXXX")

INPUT_DB="${TMP}/input_db"

if [ ! -f "${INPUT_DB}" ]; then
	mmseqs createdb \
		${INPUT_FASTA} \
		${INPUT_DB}
fi

mmseqs map \
	${INPUT_DB} \
	${TARGET_DB} \
	"${TMP}/matches" \
	${TMP} \
	--threads 32 \
	--split-memory-limit "125G"

mmseqs filterdb "${TMP}/matches" "${TMP}/top_matches" --extract-lines 1

mmseqs convertalis \
	${INPUT_DB} \
	${TARGET_DB} \
	"${TMP}/top_matches" \
	${OUTPUT_TSV} \
	--format-mode 4 \
	--format-output "query,target,qlen,tlen,fident,alnlen,mismatch,qstart,qend,tstart,tend,evalue,bits"

rm -r ${TMP}

end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"
