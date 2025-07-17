#!/bin/bash -e

## Author : Supra adapted by Abby 
## Date : 20220831

#SBATCH --job-name=AllCSO
#SBATCH --mem=16G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=19
#SBATCH --output=logs/all_cso.txt
#SBATCH --error=logs/all_cso.txt

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

#sample table containing sample name and sample index
myCSO_LibrariesFile="${current_projectPath}/${current_ver}/CSO_Samples.txt"

#file contains each antibody used in the cso and the sequence for it. There 78 in this file. 
my_Feature_Reference_file="${current_projectPath}/${current_ver}/my_FeatureReferenceFile.csv"

#####################################################################################################################################
#####################################################################################################################################

### To avoid the subshell situation, the while loop is put in a compound statement with curly brackets and an additional <read> command is used to read the column header, and there by skippng it
{
### This <read> command reads the column header i.e. it is being used to skip the 1st line in the input file which is the column header
read;
while read -r line; do


	unset -v current_CSO_SampleName current_my_LibraryFileName current_CSO_LibraryType

        current_CSO_SampleName=$(echo "${line}" | cut -f2)
        current_my_LibraryFileName=$(echo "${line}" | cut -f3)
        current_CSO_LibraryType=$(echo "${line}" | cut -f4)

	sbatch --job-name "${current_CSO_SampleName}" 045_EachSample_CSO_AntibodyCapture.sh "${current_CSO_SampleName}" "${current_projectPath}" "${current_ver}" "${current_my_LibraryFileName}" "${my_Feature_Reference_file}"

done
} < "${myCSO_LibrariesFile}"

###                                                                                                                               ###
#####################################################################################################################################
#####################################################################################################################################

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date
