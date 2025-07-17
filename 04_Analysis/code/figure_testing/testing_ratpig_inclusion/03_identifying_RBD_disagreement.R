#load the dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(stringr)
library(writexl)
library(readxl)
library(tidyseurat)

#read in seurat object
seuObj <- readRDS(file=here::here("04_Analysis", "data_objects", "04_probe", "CoVSeuratObj_VDJCSOGEX_SpecificitiesLabelled_CloneCorrected.rds"))

#make the necessary changes to the probe labels based on both the MSD data and the CITESeq data
msd <- read_xlsx(here::here("04_Analysis", "data_objects", "figure_testing", "testing_ratpig_inclusion", "CombinedMSD_Neut_Smartseq_OverlappingCITESeq.xlsx"))
mab <- read.csv(here::here("01_raw-data", "mAbMSD", "mAb_CITESeq_QC_Criteria.csv")) %>% select(!X)

mabLabelled <- mab %>%
  mutate(ProbeBinding = case_when(`WA.1` >= 1000000 & (`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto+Omi+",
                                  `WA.1` >= 1000000 & !(`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto+Omi-",
                                  !(`WA.1` >= 1000000) & (`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto-Omi+",
                                  TRUE ~ "Proto-Omi-"))

msdLabelled <- msd %>%
  group_by(clone_subgroup_id) %>%
  mutate(adj.ProtoOmi = case_when(CELL %in% mabLabelled$cell ~ mabLabelled$ProbeBinding[match(CELL, mabLabelled$cell)],
                                  CELL %in% seuObj@meta.data$CELL ~ seuObj@meta.data$adj.ProtoOmi[match(CELL, seuObj@meta.data$CELL)],
                                  TRUE ~ NA),
         Dataset = case_when(CELL %in% mabLabelled$cell ~ "CITESeq with mAb Data", TRUE ~ Dataset),
         ProbeBinding = case_when(!is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\+B\\+)|(Fully Cross-Reactive)") ~ "Proto+Omi+",
                                  !is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\-B\\+)|(BA.1-specific)") ~ "Proto-Omi+",
                                  is.na(BindingPopulation) & !is.na(adj.ProtoOmi) ~ adj.ProtoOmi,
                                  TRUE ~ "Smartseq Only")) %>%
  select(clone_subgroup_id, CELL, Dataset, ProbeBinding, BindingPopulation) %>%
  filter(ProbeBinding != "Smartseq Only") %>%
  group_by(clone_subgroup_id) %>%
  mutate(Specificity = case_when("SmartSeq" %in% Dataset ~ ProbeBinding[Dataset == "SmartSeq"][1],
                                 "CITESeq with mAb Data" %in% Dataset ~ ProbeBinding[Dataset == "CITESeq with mAb Data"][1],
                                 TRUE ~ ProbeBinding),
         Overlap = length(unique(Dataset)) > 1,
         Disagreement = length(unique(ProbeBinding)) > 1)



#reattach the labels to the CITESeq data when a clonal group is represented in the CITESeq data or the mAb data
seuObj@meta.data$adj.uncorrectedProtoOmi <- seuObj@meta.data$adj.ProtoOmi

seuObj@meta.data <- seuObj@meta.data %>%
                    mutate(adj.ProtoOmi = msdLabelled$Specificity[match(CELL, msdLabelled$Specificity)])

#save RDS
#saveRDS(seuObj,file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))
