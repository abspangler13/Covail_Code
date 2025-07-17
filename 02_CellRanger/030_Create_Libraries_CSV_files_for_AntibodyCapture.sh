#!/bin/bash -e

## Author : Supra adapted by Abby
## Date : 20223108

#SBATCH --job-name=CreateCSVFiles
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=6
#SBATCH --output=logs/make_csv.txt
#SBATCH --error=logs/make_csv.txt

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

#sample table with name of sample (containing CSOLIb at the end), the index used for library prep, and a path to a csv file for each sample names (myLib_s33...csv)
myCSO_LibrariesFile="${current_projectPath}/${current_ver}/CSO_Samples.txt"

outputLibraries_FileHeader="fastqs,sample,library_type"


#####################################################################################################################################
#####################################################################################################################################
###   Create a csv file per every CSO sample                                                                                      ###
###   i.e. for example, create file  my_libraries_S33_and_S19.csv for sample S19_316_416_Mem_WTLib & S33_316_416_Mem_CSOLib       ###

### To avoid the subshell situation, the while loop is put in a compound statement with curly brackets and an additional <read> command is used to read the column header, and there by skippng it
{
### This <read> command reads the column header i.e. it is being used to skip the 1st line in the input file which is the column header
read;
while read -r line; do


	unset -v current_CSO_SampleName current_my_LibraryFileName current_CSO_LibraryType

	current_CSO_SampleName=$(echo "${line}" | cut -f2)
	current_my_LibraryFileName=$(echo "${line}" | cut -f3)
	current_CSO_LibraryType=$(echo "${line}" | cut -f4)

	### Write the Header to the my_libraries_S1_and_S15.csv file
	echo "${outputLibraries_FileHeader}" > "${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"

	if [ "$current_CSO_SampleName" == "VRC-COVAIL-09-13-2023-CSO-Pool1" ]; then
		current_FlowCellID="AAAJCJCM5"
		raw_fastqData_folder="raw_data_other"
	else
		current_FlowCellID="HMHVTDSX7"
		raw_fastqData_folder="raw_data"
	fi

	### Append CSO Library row to the my_libraries_S1_and_S15.csv file
	echo "${current_projectPath}/${current_ver}/${raw_fastqData_folder}/${current_FlowCellID},${current_CSO_SampleName},${current_CSO_LibraryType}" >> "${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"

done
} < "${myCSO_LibrariesFile}"

#the csv file this creates just contains one line with fastq path, sample name and library type. 
###                                                                                                                               ###
#####################################################################################################################################
#####################################################################################################################################


echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date

