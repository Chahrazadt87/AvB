#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=64:mem=256gb
#PBS -e error_run_gtdbtk_de_novo_halogranum.txt
#PBS -o output_run_gtdbtk_de_novo_halogranum.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"

gtdbtk de_novo_wf \
	--genome_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/Halogranum/ \
	--out_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/GTDB-Tk_de_novo/Halogranum \
	--archaea \
	--outgroup_taxon "g__Halohasta" \
	--taxa_filter "g__Halogranum" \
	-x fna \
	--cpus 64 \
	--prefix "Halogranum" \
	--tmpdir ${TMPDIR}

