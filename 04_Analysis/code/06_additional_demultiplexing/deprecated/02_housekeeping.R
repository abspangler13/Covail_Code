#With this script, we are choosing antibodies to clone into proteins depending on our question

#load the dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(openxlsx)
library(tidyseurat)

#read in the demultiplexed object
drop.cols <- c("Population", "ProtoBA1", "ProtoOmi", "adj.Population", "TimeAdj.ProtoOmi", "SecondUMAP_1", "SecondUMAP_2","RNA.weight", "Prot.weight")

seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
df <- seuObj@meta.data %>% select(-one_of(drop.cols))
colnames(df)

# #identify a clonal threshold to choose cells from
# #Our first question: do we see appreciable maturation at the protein level of antibodies over time?
# #We expect no- there is no long-term increase in SHM of clones over time. Regardless, we want to see this by tracking
# #clones out to day 180 to see if there's any crazy stuff going on
# clones <- df[df$Infection == "N",] %>%
#   group_by(Booster, clone_subject_id, ClusterLabel, Timepoint) %>%
#   summarize(n = n()) %>%
#   ungroup() %>%
#   group_by(Booster, clone_subject_id) %>%
#   mutate(Total = sum(n),
#          TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & length(intersect(unique(Timepoint), c("Day 0", "Day 15", "Day 90"))) >= 1 ~ "Present",
#                                       TRUE ~ "Nope"),
#          Select = Total >= 15 & TimepointCorrect == "Present")
# 
# 
# day180 <- clones[clones$Select == TRUE,]
# #we will ideally choose antibodies from 1. an early timepoint and 2. day 180
# #I would also like to see antibodies from 3. a day 15 timepoint that are highly mutated. We don't expect to see
# #much maturation over time- significant increases in SHM at day 15 come from selection of more mutated B cells
# 
# #Question #2: Does infection mature B cells elicited from vaccination?
# #pull any B-cell clonal lineages from the infected group that are present both post-vax and post-infection,
# #especially people who got infected day 15-90
# clones <- df[df$Infection == "Y",] %>%
#   group_by(Booster, clone_subject_id, ClusterLabel, Timepoint) %>%
#   summarize(n = n()) %>%
#   ungroup() %>%
#   group_by(Booster, clone_subject_id) %>%
#   mutate(Total = sum(n),
#          TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & length(intersect(unique(Timepoint), c("Day 0", "Day 15", "Day 90"))) >= 1 ~ "Present",
#                                       TRUE ~ "Nope"),
#          Select = Total >= 15 & TimepointCorrect == "Present")
# 
# 
# day180Inf <- clones[clones$Select == TRUE,]

#First, we want to grab all clonal lineages present from day 0/15 and day 180 in uninfected donors
#this is less stringent than setting a 15 count threshold, which we probably won't make
clones <- df[df$Infection == "N",] %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(Booster, clone_subject_id) %>%
  mutate(Total = sum(n),
         TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & any(is.element(c("Day 0", "Day 15"),unique(Timepoint))) ~ "Present",
                                      TRUE ~ "Nope"),
         Select = TimepointCorrect == "Present")

day180 <- clones[clones$Select == TRUE,]

#Second, we want to look at donors infected days 15-90 for infected group, but we want to focus on those infected days 15-90
clones <- df[df$Infection == "Y",] %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  group_by(Booster, clone_subject_id) %>%
  mutate(Total = sum(n),
         TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & any(is.element(c("Day 0", "Day 15"),unique(Timepoint))) ~ "Present",
                                      TRUE ~ "Nope"),
         Select = TimepointCorrect == "Present")

day180Inf <- clones[clones$Select == TRUE,]

#third, just choosing all antibodies at day 180- let's label singlets or previously expanded
AllDay180 <- df %>%
  group_by(Booster, clone_subject_id) %>%
  filter(Timepoint == "Day 180") %>%
  mutate(Total = n()) %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  mutate(Filter = "Present at day 180")

#all omicron-only cells
AllOmicronOnly <- df %>% filter(adj.ProtoOmi == "Proto-Omi+") %>% mutate(Filter = "Omicron-specific cells for any group. These might be interesting to test binding for to see if they're true Omicron-specifics")

#all prototype-specific cells activated post-omicron boost
ActivatedProtoSpecific <- df %>% filter(adj.ProtoOmi == "Proto+Omi-", ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)"), Booster == "Omicron", Timepoint != "Day 0") %>% mutate(Filter = "Prototype-specific cells in an activated memory cluster for the Omicron booster. Testing these might be good for a sanity check. Note that AM3 cells are included here although it is more meaningful to focus on AM1/2 as activated cells.")

#Write an xlsx file with the antibodies we want so far
mAbList <- list()
mAbList[["UnInf_Day0or15_Day180"]] <- df[df$clone_subject_id %in% day180$clone_subject_id,] %>% group_by(clone_subject_id) %>% mutate(Total = n(), Filter = "Present at both day 0/15 and day 180")
mAbList[["Inf_Day0or15_Day180"]] <- df[df$clone_subject_id %in% day180Inf$clone_subject_id,] %>% group_by(clone_subject_id) %>% mutate(Total = n(), Filter = "Present at both day 0/15 and day 180. Filter by InfectionRange to only select those infected days 15-90")
mAbList[["AllDay180Cells"]] <- AllDay180
mAbList[["AllOmicronOnly"]] <- AllOmicronOnly
mAbList[["ActivatedProtoSpecific"]] <- ActivatedProtoSpecific

openxlsx::write.xlsx(mAbList, file= here::here("04_Analysis","data_objects","06_repertoire_analysis", "MonoclonalsToChoose_Uninfected_Infected_PresentAtDay180AndEarlier.xlsx"))
