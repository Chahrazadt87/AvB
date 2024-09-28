#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=32:mem=128gb
#PBS -e error_map_sequences_to_db.txt
#PBS -o output_map_sequences_to_db.txt

set -e

echo "The Job ID is ${PBS_JOBID}"

cd /rds/general/user/rs1521/home/

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

if [ -z "${OUTPUT_TSV}" ]; then
    echo "OUTPUT variable is not set: use qsub -v to set it"
    exit 1
fi

TMP="/rds/general/user/rs1521/ephemeral/${PBS_JOBID}"
mkdir ${TMP}

INPUT_DB="${TMP}/input_db"

if [ ! -f "${INPUT_DB}" ]; then
	mmseqs createdb \
		${INPUT_FASTA} \
		${INPUT_DB}

	# special case: map using uniprot sub db
	#mmseqs createsubdb \
	#	${INPUT_FASTA} \
	#	/rds/general/project/lms-warnecke-raw/live/UniProtKB/UniProtKB \
	#	${INPUT_DB} \
	#	--id-mode 1
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
