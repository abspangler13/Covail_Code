#!/bin/bash -e

## Author : SRG adapted by Abby
## Date : 202120831

#SBATCH --job-name=Each_GEX_CSO
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/CSO_and_WT.txt
#SBATCH --error=logs/CSO_and_WT.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load cellranger/7.2.0-dntehee

my_run_id=$1
current_projectPath=$2
current_ver=$3

current_my_LibraryFileName=$4

my_Feature_Reference_file=$5

refdata="/data/vrc_vip/Rory/Refs/refdata-gex-GRCh38-2024-A"

cr_mem=160
cr_cores=16

my_Library_file="${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"

cellranger count --id=${my_run_id} \
--libraries="${my_Library_file}" \
--feature-ref=${my_Feature_Reference_file} \
--transcriptome="${refdata}" \
--localmem=${cr_mem} \
--localcores=${cr_cores}

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date
