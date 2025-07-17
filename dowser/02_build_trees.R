suppressPackageStartupMessages(library(airr))
suppressPackageStartupMessages(library(alakazam))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dowser))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scoper))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(shazam))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(vroom))

# Bioconductor package
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(sessioninfo))

bcr_data <- readRDS(file = here::here("dowser","all_bcr_data_meta.rds"))

## resolveLight chains
bcr_data$v_gene_call = alakazam::getGene(bcr_data$v_call)
bcr_data$j_gene_call = alakazam::getGene(bcr_data$j_call)
results = resolveLightChains(bcr_data, j_call = "j_gene_call", v_call = "v_gene_call")
saveRDS(results, file = here::here("dowser","resolveLightChains.rds"))
#results <- readRDS(file = here::here("analysis","data_objects","14_dowser_redo","resolveLightChains.rds"))

# read in IMGT files in the Docker container
references <- readIMGT(dir = "/usr/local/share/germlines/imgt/human/vdj")
# reconstruct germlines
results <- createGermlines(results, references, fields = "Subject", nproc = 1)
saveRDS(results, file = here::here("dowser","createGermlines.rds"))
# results <- readRDS(file = here::here("dowser","createGermlines.rds"))

# Find all clone_ids that contain a Timepoint.num = 180
clone_ids_with_timepoint_180 <- results %>%
    filter(Timepoint.num == 180) %>%
    pull(clone_id) %>%
    unique()

# Filter results to only include those clone_ids
results <- results %>%
    filter(clone_id %in% clone_ids_with_timepoint_180)

## rename clusters
results$ClusterLabel.AS <- gsub("[\\s()]", "", results$ClusterLabel)
results$ClusterLabel.AS <- gsub(" ", "_", results$ClusterLabel.AS)

## formatClones, meta data gets removed here
my.clones.t <- formatClones(results, traits = c("Timepoint.num"), chain="HL", text_fields = c("Timepoint","Booster","adj.ProtoOmi","Infection"),
                       nproc=1, collapse = FALSE, 
                       split_light = TRUE, minseq = 3,filterstop=FALSE ) #can set split light = FALSE

saveRDS(my.clones.t, file = here::here("dowser","HL_formatClones_Timepoint.rds"))
# my.clones <- readRDS(file = here::here("dowser","HL_formatClones_Timepoint.rds"))

my.clones.s <- formatClones(results, traits = c("ProtoOmi"), chain="HL", text_fields = c("Booster","adj.ProtoOmi","Infection"),
                       nproc=1, collapse = FALSE, 
                       split_light = TRUE, minseq = 3,filterstop=FALSE ) #can set split light = FALSE

saveRDS(my.clones.s, file = here::here("dowser","HL_formatClones_Specificity.rds"))

my.clones.c <- formatClones(results, traits = c("ClusterLabel.AS"), chain="HL", text_fields = c("Booster","adj.ProtoOmi","Infection"),
                       nproc=1, collapse = FALSE, 
                       split_light = TRUE, minseq = 3,filterstop=FALSE ) #can set split light = FALSE

saveRDS(my.clones.c, file = here::here("dowser","HL_formatClones_Cluster.rds"))

## getTrees HL for Timepoint and Booster for evoltuion analysis
my.trees.t = getTrees(my.clones.t, build="raxml", 
    exec="/usr/local/bin/raxml-ng", nproc=1, partition="scaled")
saveRDS(my.trees.t, file = here::here("dowser","HL_trees_Timepoint_raxml.rds"))

# # getTrees
my.trees.s = getTrees(my.clones.s, build="raxml", 
    exec="/usr/local/bin/raxml-ng", nproc=1, partition="scaled")
saveRDS(my.trees.s, file = here::here("dowser","HL_trees_Specificity_raxml.rds"))

# getTrees
my.trees.c = getTrees(my.clones.c, build="raxml", 
    exec="/usr/local/bin/raxml-ng", nproc=1, partition="scaled")
saveRDS(my.trees.c, file = here::here("dowser","HL_trees_Cluster_raxml.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()