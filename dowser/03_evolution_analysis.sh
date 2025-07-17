#!/bin/bash -e

#SBATCH --job-name=evolution_analysis
#SBATCH --time=6:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=spanglerab@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=./dowser/logs/03_evolution_analysis.txt
#SBATCH --error=./dowser/logs/03_evolution_analysis.txt


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
    Rscript ./dowser/03_evolution_analysis.R

echo "**** Job ends ****"
date
