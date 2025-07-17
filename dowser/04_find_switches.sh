#!/bin/bash -e

#SBATCH --job-name=findSwitches
#SBATCH --mem=10G
#SBATCH --mail-user=spanglerab@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=./dowser/logs/04_find_switches.txt
#SBATCH --error=./dowser/logs/04_find_switches.txt


echo "**** Job starts ****"
date

echo "**** Skyline info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

apptainer exec \
    /data/vrc_bsc/containers/immcantation/immcantation_suite-4.5.0.sif \
    Rscript ./dowser/04_find_switches.R

echo "**** Job ends ****"
date

