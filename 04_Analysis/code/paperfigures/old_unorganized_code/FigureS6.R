#load the dependencies
library(dplyr)
library(Seurat)
library(tidyseurat)
library(stringr)
library(readxl)

set.seed(1)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObj %>% filter(ClusterLabel != "Naive")

df <- seuObj@meta.data

df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"),
         Booster = str_replace(Booster, "Omicron", "BA.1"),
         OfficialBooster = factor(OfficialBooster, levels = c("Prototype mRNA", "Prototype + Omicron BA.1 mRNA", "Omicron BA.1 mRNA")),
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")),
         InfectionRange = case_when(is.na(InfectionRange) ~ "Uninfected",
                                    !is.na(InfectionRange) ~ str_replace(InfectionRange, "Between", "Infected")))

#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#2AB673", 
               "Prototype/BA.1 mRNA" = "#1D75BC",
               "Prototype mRNA" = "#FBB042")

immunogenColors <- c("Prototype" = "#FBB042",
                     "BA.1 And Prototype" = "#1D75BC",
                     "BA.1" = "#2AB673")

rangeColors <- c("Infected Days 15-90" = "#0063B2FF",
                 "Infected Days 90-180" = "#9CC3D5FF",
                 "Uninfected" = "#ECD99f")

shortColors <- c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 #"Acute Activated" = "#F46D43",
                 "Acute Activated" = "#f08665",
                 "Intermediate" = "#E6F598",
                 "Resting IgG" = "limegreen",
                 "Resting IgA" = "#3288BD",
                 "Plasmablast-like" = "#6f2da8",
                 "Naive" = "white")

#####read in flow data
#load in the flow data and make the result
flowRaw <- read_xlsx(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Filtered_COVAILDataset_Infected.xlsx"))

flow <- flowRaw %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         Immunogen = str_replace_all(Immunogen, " \\+ ", "/"),
         Booster = str_replace_all(Booster, " \\+ ", "/")) %>%
  group_by(`Subject ID`) %>%
  filter(length(unique(Timepoint)) == 4)

flow$Timepoint <- factor(flow$Timepoint, levels = c("1", "15","90","180"))

#####
#read in evolution data from Abby
evolving <- read.csv(here::here("01_raw-data", "Evo_dat_Timepoint_uniform.csv"))
#####

#####
#define b cell clones present at different times in the response
crossDF <- df %>% filter(Infection == "Y") %>% filter(adj.ProtoOmi == "Proto+Omi+")

nonSinglets <- unique(crossDF$clone_subject_id[duplicated(crossDF$clone_subject_id) | duplicated(crossDF$clone_subject_id, fromLast=T)])
crossDF$CloneStatus <- ifelse(crossDF$clone_subject_id %in% nonSinglets, crossDF$clone_subject_id, "Singlet")

calcs <- crossDF %>%
  group_by(CloneStatus, InfectionRange, Timepoint) %>%
  summarize(n = n()) %>%
  group_by(CloneStatus, InfectionRange) %>%
  mutate(relative = case_when(unique(InfectionRange) == "Infected Days 15-90" & Timepoint %in% c("Day 90", "Day 180") ~ "Post-infection",
                              unique(InfectionRange) == "Infected Days 90-180" & Timepoint %in% c("Day 180") ~ "Post-infection",
                              TRUE ~ "Pre-infection")) %>%
  mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) & "Post-infection" %in% unique(relative) ~ "Expanded Pre-Vax, Post-Infection",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                         length(unique(Timepoint)) > 1 & length(unique(relative)) == 2 ~ "Expanded Pre- and Post-Infection",
                         length(unique(Timepoint)) > 1 & length(unique(relative)) == 1 & "Pre-infection" %in% unique(relative) ~ "Expanded Pre-Infection",
                         length(unique(relative)) == 1 & "Post-infection" %in% unique(relative) ~ "Expanded Post-Infection",
                         length(unique(relative)) == 1 & !("Post-infection" %in% unique(relative)) ~ "Single Timepoint Pre-Infection",
                         TRUE ~ "CHECK"))

crossDF$CloneStatusRefined <- calcs$lab[match(crossDF$CloneStatus, calcs$CloneStatus)]
crossDF$CloneStatusRefined <- factor(crossDF$CloneStatusRefined, levels = c("Expanded Pre-Vax, Post-Infection",
                                                                            "Expanded Pre- and Post-Infection",
                                                                            "Expanded Post-Infection",
                                                                            "Expanded Pre-Vax",
                                                                            "Expanded Pre-Infection",
                                                                            "Single Timepoint Pre-Infection",
                                                                            "Singlet"))

######do we see differences in vh usage between post-infection cells vs pre-vax
subset <- crossDF %>% filter(InfectionRange == "Infected Days 15-90", CloneStatusRefined %in% c("Expanded Post-Infection", "Single Timepoint Pre-Infection",
                                                                                                "Expanded Pre- and Post-Infection"))

subset <- subset %>% mutate(CloneStatus = case_when(CloneStatusRefined == "Expanded Post-Infection" ~ "Post-Infection",
                                                    TRUE ~ "Post-Vax"))

vh <- subset %>%
  group_by(Subject, CloneStatus, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(Subject, CloneStatus, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "CloneStatus", values_from = "Proportion")

#do stats and compare VH genes
uniqueVH <- unique(vh$v_call)

pvals <- c()
vhGene <- c()
for(i in uniqueVH){
  vhGene <- append(vhGene, i)
  vhFiltered <- vh %>% filter(v_call == i)
  
  pvals <- append(pvals, wilcox.test(vhFiltered$`Post-Infection`, vhFiltered$`Post-Vax`, paired= TRUE)$p.value)
}

pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))

#do chi-squared?
vhChi <- subset %>%
  group_by(CloneStatus, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(CloneStatus, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "CloneStatus", values_from = "Proportion")

vhChi <- as.data.frame(t(vhChi))
names(vhChi) <- unlist(vhChi[1,])
vhChi <- as.numeric(t(vhChi[-1,]))
chisq.test(vhChi)




####################################
#######also do this for activated cells after both conditions
activatedDF <- crossDF %>% 
                mutate(ActClone = case_when(ClusterLabel %in% c("Acute Activated", "Intermediate") & Timepoint == "Day 15" ~ "VaxAct",
                                            ClusterLabel %in% c("Acute Activated", "Intermediate") & Timepoint == "Day 90" & InfectionRange == "Infected Days 15-90" ~ "InfectAct",
                                            ClusterLabel %in% c("Acute Activated", "Intermediate") & Timepoint == "Day 180" & InfectionRange == "Infected Days 90-180" ~ "InfectAct",
                                            TRUE ~ NA)) %>%
                filter(!is.na(ActClone))

vh <- activatedDF %>%
  group_by(Subject, ActClone, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(Subject, ActClone, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "ActClone", values_from = "Proportion")

#do stats and compare VH genes
uniqueVH <- unique(vh$v_call)

pvals <- c()
vhGene <- c()
for(i in uniqueVH){
  vhGene <- append(vhGene, i)
  vhFiltered <- vh %>% filter(v_call == i)
  
  pvals <- append(pvals, wilcox.test(vhFiltered$VaxAct, vhFiltered$InfectAct, paired= TRUE)$p.value)
}

pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))

#do chi-squared?
vhChi <- activatedDF %>%
  group_by(ActClone, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(ActClone, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "ActClone", values_from = "Proportion")

vhChi <- as.data.frame(t(vhChi))
names(vhChi) <- unlist(vhChi[1,])
vhChi <- as.numeric(t(vhChi[-1,]))
chisq.test(vhChi)
