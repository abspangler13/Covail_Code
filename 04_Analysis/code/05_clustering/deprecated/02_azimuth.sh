#!/bin/bash -e

#$ -N azimuth
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G
#$ -m be
#$ -M rory.malek@nih.gov
#$ -o ./04_Analysis/code/05_clustering/logs/02_azimuth.txt
#$ -e ./04_Analysis/code/05_clustering/logs/02_azimuth.txt


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module unload R/4.2.1
module load R/4.1.3

Rscript ./04_Analysis/code/05_clustering/02_azimuth.R

echo "**** Job ends ****"
date