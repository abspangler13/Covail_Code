#!/bin/bash -e

#SBATCH --job-name=dowser
#SBATCH --mem=50G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/01_dowser.txt
#SBATCH --error=logs/01_dowser.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47
module load r-dowser/2.1.0-4n2ev74

Rscript 01_dowser_test.R

module unload r-dowser/2.1.0-4n2ev74

echo "**** Job ends ****"
date
