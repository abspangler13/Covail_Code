#!/bin/bash -e

#SBATCH --job-name=dsb_normalization_all
#SBATCH --mem=100G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/01_dsb_normalization_all.txt
#SBATCH --error=logs/01_dsb_normalization_all.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 01_dsb_normalization_all.R

echo "**** Job ends ****"
date