#!/bin/bash -e

## Author : SupraGajjala adapted by Abby
## Date : 20223108

#SBATCH --job-name=RunWTandCSOSamples
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/CSO_and_WT_qsub.txt
#SBATCH --error=logs/CSO_and_WT_qsub.txt

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

#contains name of csv file for each sample. That csv file contains the path for the fastq files for each type of library. 
myCSO_LibrariesFile="${current_projectPath}/${current_ver}/CSO_Libraries.txt"

my_Feature_Reference_file="${current_projectPath}/${current_ver}/my_FeatureReferenceFile.csv"

#####################################################################################################################################
#####################################################################################################################################
###   Create a csv file per every CSO sample                                                                                      ###
###   i.e. for example, create file  my_libraries_S33_and_S19.csv for sample S19_316_416_Mem_WTLib & S33_316_416_Mem_CSOLib       ###

### To avoid the subshell situation, the while loop is put in a compound statement with curly brackets and an additional <read> command is used to read the column header, and there by skippng it
{
### This <read> command reads the column header i.e. it is being used to skip the 1st line in the input file which is the column header
read;
while read -r line; do


	unset -v current_run_id current_my_LibraryFileName

	current_run_id=$(echo "${line}" | cut -f1)
	current_my_LibraryFileName=$(echo "${line}" | cut -f6)

	sbatch --job-name "${current_run_id}" 065_EachSample_WT_and_CSO.sh "${current_run_id}" "${current_projectPath}" "${current_ver}" "${current_my_LibraryFileName}" "${my_Feature_Reference_file}"

done
} < "${myCSO_LibrariesFile}"

###                                                                                                                               ###
#####################################################################################################################################
#####################################################################################################################################

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date