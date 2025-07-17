#!/bin/bash -e

#SBATCH --job-name=vdj1
#SBATCH --mem=100G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/01_vdj_all.txt
#SBATCH --error=logs/01_vdj_all.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 01_vdj_all.R

echo "**** Job ends ****"
date