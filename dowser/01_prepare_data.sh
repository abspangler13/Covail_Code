#!/bin/bash -e

#SBATCH --job-name=prepare_data
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=spanglerab@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=./dowser/logs/01_prepare_data.txt
#SBATCH --error=./dowser/logs/01_prepare_data.txt


echo "**** Job starts ****"
date

echo "**** Skyline info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURMD_NODENAME}"
pwd

# Load the r-tidyseurat module version 0.8.0 for compatibility with the analysis script
module load r/4.3.0-a22xr47

Rscript ./dowser/01_prepare_data.R

echo "**** Job ends ****"
date
