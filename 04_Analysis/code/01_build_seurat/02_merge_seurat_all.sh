#!/bin/bash -e

#SBATCH --job-name=merge_seurat_all
#SBATCH --mem=100G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --output=logs/merge_seurat_all.txt
#SBATCH --error=logs/merge_seurat_all.txt

echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load r/4.3.0-a22xr47

Rscript 02_merge_seurat_all.R

echo "**** Job ends ****"
date