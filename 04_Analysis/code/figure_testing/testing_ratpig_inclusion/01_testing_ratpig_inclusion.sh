#!/bin/bash -e

#SBATCH --job-name=01_testing_ratpig_inclusion
#SBATCH --mem=50G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/01_testing_ratpig_inclusion.txt
#SBATCH --error=logs/01_testing_ratpig_inclusion.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 01_testing_ratpig_inclusion.R

echo "**** Job ends ****"
date
