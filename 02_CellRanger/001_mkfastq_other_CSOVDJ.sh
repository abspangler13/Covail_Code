#!/bin/bash -e

## Adapted from Supra 
## Date : 110723

#SBATCH --job-name=myMakeFastq_ExtraCSOVDJ
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/mkfastq_other_CSOVDJ.txt
#SBATCH --error=logs/mkfastq_other_CSOVDJ.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load cellranger/7.2.0-dntehee
module load bcl2fastq2/2.20.0.422-orocbiu


run_name="/data/vrc_bsc/Rory/Covail/01_raw-data/230912_VL00414_28_AAAJCJCM5"

current_projectPath="/data/vrc_bsc/Rory/Covail/"

current_ver="02_CellRanger"

sample_sheet="${current_projectPath}/${current_ver}/SampleSheet_CSOVDJ_other.csv"

output_path="${current_projectPath}/${current_ver}/raw_data_other"

cr_cores=16
cr_mem=160

# had to make changes cuz of NovaSeq upgrade
# https://kb.10xgenomics.com/hc/en-us/articles/6748964203149-Will-mkfastq-be-affected-by-the-NovaSeq-control-software-upgrade-to-v1-8-0-?source=answerbot 
cellranger mkfastq \
--run "${run_name}" \
--csv "${sample_sheet}" \
--output-dir "${output_path}" \
--rc-i2-override=true \
--localcores=${cr_cores} \
--localmem=${cr_mem}

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date