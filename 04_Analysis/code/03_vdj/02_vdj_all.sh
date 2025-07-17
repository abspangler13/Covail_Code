#!/bin/bash -e

#SBATCH --job-name=vdj2
#SBATCH --mem=100G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/02_vdj_all.txt
#SBATCH --error=logs/02_vdj_all.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 02_vdj_all.R

echo "**** Job ends ****"
date