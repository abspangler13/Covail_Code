#!/bin/bash -e

## Author : Supra Gajjala adapted by Abby
## Date : 202230831

#SBATCH --job-name=CreateCustomCSOGEX
#SBATCH --mem=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=19
#SBATCH --output=logs/create_csv_cso_custom.txt
#SBATCH --error=logs/create_csv_cso_custom.txt

module load cellranger/7.2.0-dntehee

current_projectPath="/data/vrc_bsc/Rory/Covail"

current_ver="02_CellRanger"

#text file containing sample names and names of both the WT and CSO samples names. Then it contains a filename for the my_libraries...csv for each sample. 
myCSO_LibrariesFile="${current_projectPath}/${current_ver}/CSO_Libraries.txt"

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


	unset -v current_run_id current_WT_SampleName current_CSO_SampleName current_WT_LibraryType current_CSO_LibraryType current_my_LibraryFileName

	current_run_id=$(echo "${line}" | cut -f1)
	current_WT_SampleName=$(echo "${line}" | cut -f2)
	current_CSO_SampleName=$(echo "${line}" | cut -f3)
	current_WT_LibraryType=$(echo "${line}" | cut -f4)
	current_CSO_LibraryType=$(echo "${line}" | cut -f5)
	current_my_LibraryFileName=$(echo "${line}" | cut -f6)
	

	### Write the Header to the my_libraries_S1_and_S15.csv file
	if [ "$current_run_id" == "VRC-COVAIL-09-13-2023-Pool1" ]; then
		current_FlowCellID_CSO="AAAJCJCM5"
		current_FlowCellID_GEX="AAF5GCGM5"
		raw_fastqData_folder_CSO="raw_data_other"
		raw_fastqData_folder_GEX="raw_data_other_gex"
	else
		current_FlowCellID_CSO="HMHVTDSX7"
		current_FlowCellID_GEX="HMHVTDSX7"
		raw_fastqData_folder_GEX="raw_data"
		raw_fastqData_folder_CSO="raw_data"
	fi

	### Adds headers to lib file name
	echo "${outputLibraries_FileHeader}" > "${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"

	### Append WTO Library row to the my_libraries_S1_and_S15.csv file
	echo "${current_projectPath}/${current_ver}/${raw_fastqData_folder_GEX}/${current_FlowCellID_GEX},${current_WT_SampleName},${current_WT_LibraryType}" >> "${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"

	### Append CSO Library row to the my_libraries_S1_and_S15.csv file
	echo "${current_projectPath}/${current_ver}/${raw_fastqData_folder_CSO}/${current_FlowCellID_CSO},${current_CSO_SampleName},${current_CSO_LibraryType}" >> "${current_projectPath}/${current_ver}/${current_my_LibraryFileName}"


done
} < "${myCSO_LibrariesFile}"

#creates file named liked this my_libraries_S33_and_S17.csv 
###                                                                                                                               ###
#####################################################################################################################################
#####################################################################################################################################


