#!/bin/bash -e

## Author : Supra Adapted byAbby
## Date : 20220831

#SBATCH --job-name=AllVDJ
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks=19
#SBATCH --partition=all
#SBATCH --output=logs/all_vdj.txt
#SBATCH --error=logs/all_vdj.txt

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

#file with sample name and sample index used during library prep. 16 samples total. 
myVDJSamplesFile="${current_projectPath}/${current_ver}/VDJ_Samples.txt"

myVDJ_SAMPLES=$(tail -n+2 ${myVDJSamplesFile} | cut -f2)

#loop that will iterate through each sample name/ID from the VDJ_Samples.txt file and qsub a job
for each_sampleID in ${myVDJ_SAMPLES}
do
	echo -e "${each_sampleID}"

	#qsub 015_EachSample_VDJcount.sh "${each_sampleID}" "${current_projectPath}" "${current_ver}"

	sbatch -J "${each_sampleID}" 015_EachSample_VDJcount.sh "${each_sampleID}" "${current_projectPath}" "${current_ver}"

done

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date

