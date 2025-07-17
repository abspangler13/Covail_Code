#!/bin/bash -e

#SBATCH --job-name=reclustering
#SBATCH --mem=50G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/03_reclustering_04res.txt
#SBATCH --error=logs/03_reclustering_04res.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 03_reclustering_04res.R

echo "**** Job ends ****"
date
