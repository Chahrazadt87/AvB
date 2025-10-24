#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=32:mem=256gb
#PBS -e error_run_gtdbtk_halo.txt
#PBS -o output_run_gtdbtk_halo.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"

mkdir -p ${TMPDIR}
mkdir ${TMPDIR}/scratch_dir

gtdbtk classify_wf \
	--genome_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/archaea/ \
	--out_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/GTDB-Tk_outputs/archaea/ \
	--skip_ani_screen \
	-x fna \
	--cpus 32 \
	--prefix "Halo" \
	--tmpdir ${TMPDIR} \
	--scratch_dir ${TMPDIR}/scratch_dir \
	--pplacer_cpus 1

