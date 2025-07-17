#!/bin/bash -e

#SBATCH --job-name=Immcantation
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=rory.malek@nih.gov
#SBATCH --mail-type=BEGIN,END
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/immcantation.txt
#SBATCH --error=logs/immcantation.txt
#SBATCH --array=1-30

# code from:
# Immcantation 10x vdj tutorial https://immcantation.readthedocs.io/en/stable/tutorials/10x_tutorial.html
# Immcantation igblasthttps://changeo.readthedocs.io/en/stable/examples/10x.html?highlight=LOCUS#splitting-into-separate-light-and-heavy-chain-files 
# modeled after this file: smb://vrc-fls-b.niaid.nih.gov/IMC/Abby/VDJ analysis/To do clonal analysis using changeo.docx


echo "**** Job starts ****"
date

echo "**** LOCUS info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
pwd

module load immcantation
module load r/4.3.0-a22xr47

# get correct line from sample table
SAMPLE_LINE=$(head -n $SLURM_ARRAY_TASK_ID vdj_files.txt | tail -1)
FASTA=$(echo "${SAMPLE_LINE}" | cut -f2 -d$'\t')
ANN=$(echo "${SAMPLE_LINE}" | cut -f3 -d$'\t')
SAMPLE_ID=$(echo "${SAMPLE_LINE}" | cut -f1 -d$'\t')


#to specify the immcantation container, we'll need to specifically reference the .sif file and not just the folder it's contained in. Locus previously stored this as $IMM_SING
IMMCAN_APP=/data/apps/software/containers/immcantation/immcantation_suite-4.4.0.sif

echo ${SAMPLE_LINE}

#assign VDJ genes by IGBlast
apptainer exec $IMMCAN_APP AssignGenes.py igblast \
    -s $FASTA \
    -b /usr/local/share/igblast  \
    -o "${SAMPLE_ID}_filtered_contig_igblast.fmt7" \
    --organism human \
    --loci ig \
    --format blast
# produces this file All_filtered_contig_igblast.fmt7. Alignes your sequences to germline. determines which genes the sequence comes from. 

#please note that the x in "10x" should be lowercase!!
#make database
apptainer exec $IMMCAN_APP MakeDb.py igblast \
    -i "${SAMPLE_ID}_filtered_contig_igblast.fmt7" \
    -s $FASTA \
    -r /usr/local/share/germlines/imgt/human/vdj \
    --10x ${ANN}
# produces this file filtered_contig_igblast_db-pass.tsv. Takes infor from AssignGenes and puts it in a readable format

#parse database
apptainer exec $IMMCAN_APP /usr/local/bin/ParseDb.py select \
    -d "${SAMPLE_ID}_filtered_contig_igblast_db-pass.tsv" \
    -f locus \
    -u "IGH" \
    --logic all \
    --regex \
    --outname "${SAMPLE_ID}_heavy"
# produces this file S8_heavy_parse-select.tsv. This pulls out all the sequences that aligned to heavy chain

#set clonal labels based on CDRH3 distance + junction length, etc
apptainer exec $IMMCAN_APP /usr/local/bin/DefineClones.py \
    -d "${SAMPLE_ID}_heavy_parse-select.tsv"  \
    --outname "${SAMPLE_ID}_heavy" \
    --act set \
    --model ham \
    --norm len \
    --dist 0.15 #parameter for determining at what point we consider things clones of eachother
# produces S8_clone-pass.tsv.  Takes heavy chains determines which ones are clonaly related to each other. only based on heavy chains. determines they're clones despite hypermutation

#do for light chain now
apptainer exec $IMMCAN_APP /usr/local/bin/ParseDb.py select \
    -d "${SAMPLE_ID}_filtered_contig_igblast_db-pass.tsv" \
    -f locus \
    -u "IG[LK]" \
    --logic all \
    --regex \
    --outname "${SAMPLE_ID}_light"
# produces S8_light_parse-select.tsv. Pull out all the light chains

#from here https://changeo.readthedocs.io/en/stable/examples/10x.html?highlight=light_cluster
apptainer exec $IMMCAN_APP /usr/local/bin/light_cluster.py \
    -d "${SAMPLE_ID}_heavy_clone-pass.tsv" \
    -e "${SAMPLE_ID}_light_parse-select.tsv" \
    -o "${SAMPLE_ID}_heavy_light_clone-pass.tsv"
#produces S8_light_clone-pass.tab Redefine clones based on light chain info. Should call this HC_LC_clone-pass

#append germline sequences to dataset
apptainer exec $IMMCAN_APP /usr/local/bin/CreateGermlines.py \
    -d "${SAMPLE_ID}_heavy_light_clone-pass.tsv" \
    -g dmask \
    --cloned \
    -r /usr/local/share/germlines/imgt/human/vdj

echo "**** Job ends ****"
echo "Job id: ${JOB_ID}"
date