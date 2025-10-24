#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8:mem=128gb
#PBS -e error_run_gtdbtk_pontibacillus.txt
#PBS -o output_run_gtdbtk_pontibacillus.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"

gtdbtk classify_wf \
	--genome_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/bacteria/ \
	--out_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/GTDB-Tk_outputs/bacteria/ \
	--skip_ani_screen \
	-x fna \
	--cpus 8 \
	--prefix "Pontibacillus" \
	--tmpdir ${TMPDIR}

