#!/bin/bash -e

#SBATCH --job-name=addingSHM
#SBATCH --mem=20G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/01a_adding_SHM.txt
#SBATCH --error=logs/01a_adding_SHM.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47
module load r-shazam/1.2.0-r2m4hmu

Rscript 01a_adding_SHM.R

module unload r-shazam/1.2.0-r2m4hmu

echo "**** Job ends ****"
date