#!/bin/bash -e

## Author : Supra Adated by Abby
## Date : 20220831

#SBATCH --job-name=EachVDJ
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/each_vdj.txt
#SBATCH --error=logs/each_vdj.txt

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

if [ "$my_run_id" == "VRC-COVAIL-09-13-2023-VDJ-Pool1" ]; then
	fastqs_path="/data/vrc_bsc/Rory/Covail/02_CellRanger/raw_data_other"
else
	fastqs_path="/data/vrc_bsc/Rory/Covail/02_CellRanger/raw_data"
fi

###refdata="/hpcdata/vrc/vrc1_data/douek_lab/reference_sets/tenX/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0"
refdata="/data/vrc_bsc/refdata_9WWCAVD/tenX/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0"

cr_mem=160
cr_cores=16

if [ "$my_run_id" == "VRC-COVAIL-09-13-2023-VDJ-Pool1" ]; then
	params="AAAJCJCM5"
else
	params="HMHVTDSX7"
fi

echo "Running job ${my_run_id} from flow cell ${params}. FASTQ files are found at this path: ${fastqs_path}."

cellranger vdj --id=${my_run_id} \
--reference="${refdata}" \
--fastqs="${fastqs_path}" \
--sample="${sample_name}" \
--localcores=${cr_cores} \
--localmem=${cr_mem} \
--project="${params}"

#produces folder for each sample that starts with sample name then contains VDJLib. Example: S1_A1_408_Mem_Mem_B_cells_VDJLib

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date


