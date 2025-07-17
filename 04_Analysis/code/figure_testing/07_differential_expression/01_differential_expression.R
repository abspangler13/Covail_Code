library(Seurat)
library(dplyr)
library(tidyseurat)
library(ggplot2)
library(ggalluvial)
library(here)

#load in the data
seuObjNaive <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObjNaive %>% filter(ClusterLabel != "Naive")
df <- seuObj@meta.data #creates a dataframe of the metadata

df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"))


df$ClusterLabel <- factor(df$ClusterLabel, levels = c("AM1 (Activated)",
                                                      "AM2 (Intermediate)",
                                                      "AM3 (Atypical)",
                                                      "Activated IgA Memory",
                                                      "Resting IgA Memory",
                                                      "Resting IgG Memory",
                                                      "Resting IgG Memory 2",
                                                      "Unclear"))

# df$ClusterLabel <- factor(df$ClusterLabel, levels = c("C1-Activated",
#                                                       "C2-Atypical",
#                                                       "C4-Activated IgA Memory",
#                                                       "C3-Intermediate",
#                                                       "C6-Resting IgA Memory",
#                                                       "C8-Resting IgG Memory",
#                                                       "C9-Resting IgG Memory 2",
#                                                       "Unclear"))
# 
# df$OfficialInfection <- ifelse(df$Infection == "Y", "Infected", "Uninfected")
# df$OfficialInfection <- factor(df$OfficialInfection, levels = c("Uninfected",
#                                                                 "Infected"))
# 
# df$InfectionRange <- ifelse(is.na(df$InfectionRange), "Uninfected", df$InfectionRange)

#set colors
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#045275")

phenoColors <- c("AM1 (Activated)" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 "AM2 (Intermediate)" = "#F46D43",
                 "AM3 (Atypical)" = "#FDAE61",
                 "Activated IgA Memory" = "#FEE08B",
                 "Resting IgA Memory" = "#E6F598",
                 "Resting IgG Memory" = "#ABDDA4",
                 "Resting IgG Memory 2" = "#66C2A5",
                 "Unclear" = "#3288BD",
                 "Naive" = "#6f2da8")

########
#Run DE analysis
########
check <- NormalizeData(seuObjNaive)
Idents(check) <- "ClusterLabel"
DefaultAssay(check) <- "RNA"

#differences between resting compartments
markers <- FindMarkers(check, ident.1 = c("Resting IgA Memory", "Activated IgA Memory"), ident.2 = c("Resting IgG Memory", "Resting IgG Memory 2")) %>%
            filter(abs(avg_log2FC) > 1.5, !grepl("IGH", rownames(.)))

