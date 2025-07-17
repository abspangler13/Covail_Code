#!/bin/bash -e

## Author : Supra adapted by Abby
## Date : 20220831

#SBATCH --job-name=AllGEX_Count
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks=19
#SBATCH --partition=all
#SBATCH --output=logs/all_samples_GEX_count.txt
#SBATCH --error=logs/all_samples_GEX_count.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load cellranger/7.2.0-dntehee

current_projectPath="/data/vrc_bsc/Rory/Covail"

current_ver="02_CellRanger"

#sample table that contains sample name and index used during library prep. Sample names contain "WT"
myCurrentSamplesFile="${current_projectPath}/${current_ver}/WT_Samples.txt"

my_SAMPLES=$(tail -n+2 ${myCurrentSamplesFile} | cut -f2)

#loop to submit job for each sample
for each_sampleID in ${my_SAMPLES}
do
	echo -e "${each_sampleID}"

	#qsub 015_EachSample_VDJcount.sh "${each_sampleID}" "${current_projectPath}" "${current_ver}"

	sbatch -J "${each_sampleID}" 025_EachSample_WTcount.sh "${each_sampleID}" "${current_projectPath}" "${current_ver}"

done


echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date
