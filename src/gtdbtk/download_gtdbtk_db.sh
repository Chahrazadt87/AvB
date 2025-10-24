#!/bin/bash
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=1:mem=8gb
#PBS -e error_download_gtdbtk_db.txt
#PBS -o output_download_gtdbtk_db.txt

cd /rds/general/user/rs1521/home/

. load_conda.sh
conda activate gtdbtk-2.4.0

mkdir -p /rds/general/project/lms-warnecke-raw/live/gtdbtk-2.4.0

download-db.sh /rds/general/project/lms-warnecke-raw/live/gtdbtk-2.4.0

