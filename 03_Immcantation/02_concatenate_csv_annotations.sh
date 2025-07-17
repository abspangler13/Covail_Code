#!/bin/bash -e

#SBATCH --job-name=make_csv_annotations
#SBATCH --mem=20G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/make_ann.txt
#SBATCH --error=logs/make_ann.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 02_concatenate_csv_annotations.R

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date