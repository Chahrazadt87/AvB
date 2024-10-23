#!/bin/bash
#PBS -l walltime=01:30:00
#PBS -l select=1:ncpus=16:mem=512gb
#PBS -e error_run_foldseek_test.txt
#PBS -o output_run_foldseek_test.txt

# Running on Imperial's HPC.

set -e

start=$(date +%s)

BASE=/rds/general/user/rs1521/home
TMP=/rds/general/ephemeral/user/rs1521/ephemeral
AFDB=/rds/general/project/lms-warnecke-raw/live/afdb50

WORKDIR=$TMP/run_foldseek
mkdir -p $WORKDIR

cd $BASE

# Using a more recent foldseek version than the one available on conda
export PATH=/rds/general/user/rs1521/home/foldseek_ef4e960/build/bin/:$PATH

foldseek easy-search \
	$BASE/chahrazad/hydrolase_2090_AF3.cif \
	$AFDB/afdb50 \
	$WORKDIR/AF-J2ZZK6-F1-model_v4_results.tsv \
	$TMP \
	--threads 16 \
	--prefilter-mode 1 \
	--cluster-search 1 \
       	--format-mode 4 \
        --format-output query,target,taxid,taxname,taxlineage,fident,alnlen,tstart,tend,tlen,evalue,bits,qtmscore,rmsd

cp \
	$WORKDIR/AF-J2ZZK6-F1-model_v4_results.tsv \
	$BASE/chahrazad/hydrolase_2090_AF3_foldseek_results.tsv

end=$(date +%s)
echo "Elapsed Time: $(($end-$start)) seconds"
