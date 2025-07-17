#demultiplexing infected cohort and labelling phenotypic classifications from Sarah's paper

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

#load the data- use resolution 04
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "05_clustering", "res_04", "COVAIL_ReclusteredAzimuth_04res.rds"))

#let's first demultiplex by infection
seuObj@meta.data <- seuObj@meta.data %>%
                    mutate(InfectionTimepoint = case_when(Infection == "N" ~ "Uninfected",
                                                 Timepoint %in% c("Day 90", "Day 180") & Subject %in% c(5755544848,
                                                                                                        5750564848,
                                                                                                        5657574848,
                                                                                                        5556484948,
                                                                                                        5553554848,
                                                                                                        5356534848,
                                                                                                        5249544848,
                                                                                                        5050484948,
                                                                                                        5048544848) ~ "Post-Infection",
                                                 Timepoint == "Day 180" & Subject %in% c(5653544848, 
                                                                                         5456544848, 
                                                                                         5054574848, 
                                                                                         4950544848,
                                                                                         4856554848,
                                                                                         5351564848) ~ "Post-Infection",
                                                 TRUE ~ "Pre-Infection"),
                                 InfectionRange = case_when(Infection == "N" ~ NA,
                                                            Subject %in% c(5755544848,
                                                                           5750564848,
                                                                           5657574848,
                                                                           5556484948,
                                                                           5553554848,
                                                                           5356534848,
                                                                           5249544848,
                                                                           5050484948,
                                                                           5048544848) ~ "Between Days 15-90",
                                                            Subject %in% c(5653544848, 
                                                                           5456544848, 
                                                                           5054574848, 
                                                                           4950544848,
                                                                           4856554848) ~ "Between Days 90-180",
                                                            TRUE ~ "Check again"))

#bring in a metadata file and finish up the rest of the unblinding
unblinding <- read.csv(here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "metadata", "VRC_infections_Details_01MAR2024_Stg1_2_3.csv")) %>%
              select(Subject_ID, PANGOLIN_SWAB1, TRUNCLIN_SWAB1, datinf, source,ENRLDATE) %>%
              rename(Subject=Subject_ID) %>%
              mutate(Subject = as.character(Subject))

#we will now need to ascribe labels to the different populations (and combine some):
seuObj@meta.data <- seuObj@meta.data %>%
                      mutate(ClusterLabel = case_when(seurat_clusters %in% c(0,3,8) ~ "Resting IgG",
                                                      seurat_clusters == 1  ~ "Intermediate",
                                                      seurat_clusters == 2 ~ "Acute Activated",
                                                      seurat_clusters == 4 ~ "Resting IgA",
                                                      seurat_clusters == 5 ~ "Atypical",
                                                      seurat_clusters == 6 ~ "Intermediate",
                                                      seurat_clusters == 7 ~ "Naive",
                                                      seurat_clusters == 9 ~ "Plasmablast-like",
                                                      TRUE ~ "Check code pls"),
                             EnrollmentDate = unblinding$ENRLDATE[match(Subject, unblinding$Subject)],
                             InfectionDate = unblinding$datinf[match(Subject, unblinding$Subject)],
                             InfectionLineage = unblinding$PANGOLIN_SWAB1[match(Subject, unblinding$Subject)],
                             InfectionTrunc = unblinding$TRUNCLIN_SWAB1[match(Subject, unblinding$Subject)],
                             InfectionIndentMethod = unblinding$source[match(Subject, unblinding$Subject)])

seuObj@meta.data$ClusterLabel <- factor(seuObj@meta.data$ClusterLabel, levels = c("Acute Activated",
                                                      "Intermediate",
                                                      "Atypical",
                                                      "Resting IgG",
                                                      "Resting IgA",
                                                      "Plasmablast-like",
                                                      "Naive"))

#save RDS
saveRDS(seuObj,file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))
