#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l select=1:ncpus=64:mem=256gb
#PBS -e error_run_gtdbtk_de_novo_strep.txt
#PBS -o output_run_gtdbtk_de_novo_strep.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"

gtdbtk de_novo_wf \
	--genome_dir /rds/general/user/rs1521/home/Streptococcus/fna_files/ \
	--out_dir /rds/general/user/rs1521/home/Streptococcus/GTDB-Tk_de_novo/ \
	--bacteria \
	--outgroup_taxon "g__Lactococcus" \
	--taxa_filter "g__Streptococcus" \
	-x gz \
	--cpus 64 \
	--prefix "Streptococcus" \
	--tmpdir ${TMPDIR}

