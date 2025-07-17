#!/bin/bash -e

#SBATCH --job-name=completingAnnotations
#SBATCH --mem=20G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/02_completingannotations.txt
#SBATCH --error=logs/02_completingannotations.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 02_completing_annotations.R

echo "**** Job ends ****"
date