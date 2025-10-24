#!/bin/bash
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=64:mem=256gb
#PBS -e error_run_gtdbtk_de_novo_ponti.txt
#PBS -o output_run_gtdbtk_de_novo_ponti.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}"

gtdbtk de_novo_wf \
	--genome_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/Pontibacillus/ \
	--out_dir /rds/general/user/rs1521/home/chahrazad/Atanasova_genomes/GTDB-Tk_de_novo/Pontibacillus/ \
	--bacteria \
	--outgroup_taxon "f__Halobacillaceae" \
	--taxa_filter "g__Pontibacillus" \
	-x fna \
	--cpus 64 \
	--prefix "Pontibacillus" \
	--tmpdir ${TMPDIR}

