library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(stringr)
library(writexl)
library(readxl)
library(tidyseurat)

#########write new seurat object and export metadata for Flavio to use
seuObj <- readRDS(file=here::here("04_Analysis", "data_objects", "04_probe", "CoVSeuratObj_VDJCSOGEX_SpecificitiesLabelled_CloneCorrected.rds"))

#########make the necessary changes to the probe labels based on both the MSD data and the CITESeq data
msd <- read_xlsx(here::here("01_raw-data", "CombinedMSD_Neut_Smartseq_OverlappingCITESeq.xlsx"))
mab <- read.csv(here::here("01_raw-data", "mAbMSD", "mAb_CITESeq_QC_Criteria.csv")) %>% select(!X)

#defines positivity for mabs based upon different ECL unit thresholds
mabLabelled <- mab %>%
  mutate(ProbeBinding = case_when(`WA.1` >= 1000000 & (`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto+Omi+",
                                `WA.1` >= 1000000 & !(`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto+Omi-",
                                !(`WA.1` >= 1000000) & (`BA.1` >= 500000 | `XBB` >= 500000) ~ "Proto-Omi+",
                                TRUE ~ "Proto-Omi-"))

#start the correction
msdLabelled <- msd %>%
  group_by(clone_subgroup_id) %>%
  mutate(adj.ProtoOmi = case_when(CELL %in% mabLabelled$cell ~ mabLabelled$ProbeBinding[match(CELL, mabLabelled$cell)],
                                  CELL %in% seuObj@meta.data$CELL ~ seuObj@meta.data$adj.ProtoOmi[match(CELL, seuObj@meta.data$CELL)],
                                  TRUE ~ adj.ProtoOmi), #when present in mAb data, set CITESeq call to match the mAbs
         Dataset = case_when(CELL %in% mabLabelled$cell ~ "CITESeq with mAb Data", TRUE ~ Dataset), #label as CITESeq cells we have mAb data for
         ProbeBinding = case_when(!is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\+B\\+)|(Fully Cross-Reactive)") ~ "Proto+Omi+",
                                  !is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\-B\\+)|(BA.1-specific)") ~ "Proto-Omi+",
                                  is.na(BindingPopulation) & !is.na(adj.ProtoOmi) ~ adj.ProtoOmi,
                                  TRUE ~ "Smartseq Only")) %>% #set a universal probe-binding variable to compare ratpig/citeseq/mab calls
  select(clone_subgroup_id, CELL, Dataset, ProbeBinding, BindingPopulation) %>%
  filter(ProbeBinding != "Smartseq Only") %>% #remove cells only represented in our smartseq data and not our citeseq data
  group_by(clone_subgroup_id) %>%
  mutate(Specificity = case_when("SmartSeq" %in% Dataset ~ ProbeBinding[Dataset == "SmartSeq"][1],
                                 "CITESeq with mAb Data" %in% Dataset ~ ProbeBinding[Dataset == "CITESeq with mAb Data"][1],
                                 TRUE ~ ProbeBinding)) #default to the binding label being the one we have actual binding data for, otherwise carry on :)

#add the changes back into the seurat object
seuObj@meta.data$uncorrected.adj.ProtoOmi <- seuObj@meta.data$adj.ProtoOmi
seuObj@meta.data$adj.ProtoOmi <- msdLabelled$Specificity[match(seuObj@meta.data$CELL, msdLabelled$CELL)]

#test for proper clonal correction
test <- seuObj@meta.data %>% group_by(clone_subject_id) %>% mutate(Check = length(unique(adj.ProtoOmi)) > 1) %>% select(clone_subject_id, adj.ProtoOmi, uncorrected.adj.ProtoOmi, Check)

print("Presence of non-matching clonal specificity labels:")
print(nrow(test[test$Check,]))

#save the corrected seurat object
saveRDS(seuObj, here::here("04_Analysis", "data_objects", "04_probe", "BindingDataCorrected_COVAIL_SeuratObj.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

