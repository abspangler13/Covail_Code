#!/bin/bash -e

## Author : Supra adapted by Abby
## Date : 20220831

#SBATCH --job-name=EachSampleGEX_count
#SBATCH --mem=16G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=6
#SBATCH --output=logs/each_GEX_count.txt
#SBATCH --error=logs/each_GEX_count.txt

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
sample_name=$1

current_projectPath=$2
current_ver=$3

if [ "$my_run_id" == "VRC-COVAIL-09-13-2023-GEX-Pool1" ]; then
	fastqs_path="/data/vrc_bsc/Rory/Covail/02_CellRanger/raw_data_other_gex"
else
	fastqs_path="/data/vrc_bsc/Rory/Covail/02_CellRanger/raw_data"
fi

refdata="/data/vrc_bsc/Rory/Refs/refdata-gex-GRCh38-2024-A"

if [ "$my_run_id" == "VRC-COVAIL-09-13-2023-GEX-Pool1" ]; then
	params="AAF5GCGM5"
else
	params="HMHVTDSX7"
fi

echo "Running ${my_run_id}. Fastq path is ${fastqs_path}, flow cell ${params}."

cr_mem=60
cr_cores=6

cellranger count --id=${my_run_id} \
--transcriptome="${refdata}" \
--fastqs="${fastqs_path}" \
--sample="${sample_name}" \
--chemistry=fiveprime \
--localcores=${cr_cores} \
--localmem=${cr_mem} \
--project="${params}"


echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date