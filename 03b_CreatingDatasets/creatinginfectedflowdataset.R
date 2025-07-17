library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

#load in the flow data run on XBB
metadata <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx"))
metadata2 <- read.csv(here::here("01_raw-data", "FlowData","bcell_unblinding.csv"))

flow <- read_xlsx(here::here("01_raw-data", "FlowData","241004_CombinedXBBAndDeltaDatasets.xlsx"))

#identify infected/oos boosted people
oosBoost <- unique(flow$`Subject ID`[flow$oosboost_flag == "Y" & !is.na(flow$infect_flag)])
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

#continue editing data
flow2 <- flow %>%
  filter(!is.na(`Specimen ID`)) %>%
  mutate(`Time point Guess` = as.character(`Time point Guess`),
         Timepoint = case_when(timepoint == "Day 1" ~ "1",
                               timepoint == "Day 15" ~ "15",
                               timepoint == "Day 57" ~ "57",
                               timepoint == "Day 71" ~ "71",
                               timepoint == "Day 91" ~ "90",
                               timepoint == "Day 181" ~ "180",
                               TRUE ~ `Time point Guess`),
         CompleteBoost = case_when(Dataset == "XBB Panel" ~ metadata$Treatment[match(`Subject ID`, metadata$`Subject ID`)],
                                   Dataset == "Delta Panel" ~ metadata2$treatment[match(`Specimen ID`, metadata2$sn)]),
         Company = str_extract(CompleteBoost, "(Sanofi)|(Pfizer)|(Moderna)"),
         Treatment = ifelse(treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", treatment),
         `Vaccine Platform` = case_when(Dataset == "XBB Panel" ~ "Not Applicable",
                                        `Vaccine Platform` == "Sanofi" ~ "Protein",
                                        `Vaccine Platform` == "Pfizer" ~ "mRNA",
                                        `Vaccine Platform` == "Moderna" ~ "mRNA"),
         Treatment = case_when(Dataset == "XBB Panel" ~ Treatment,
                               Dataset == "Delta Panel" ~ paste0(Treatment, " ", `Vaccine Platform`),
                               TRUE ~ "ERROR"),
         Booster = case_when(Treatment %in% c("1 Dose Prototype (Moderna)", "Wildtype/Prototype (Pfizer 1)", "Prototype mRNA") ~ "Prototype mRNA",
                             Treatment %in% c("1 Dose Omicron (Moderna)", "Omicron (Pfizer 1)", "Omicron mRNA") ~ "Omicron BA.1 mRNA",
                             Treatment %in% c("1 Dose Omicron + Prototype (Moderna)", "Omicron + Prototype mRNA", "Omicron + Wildtype/Prototype (Pfizer 1)") ~ "Prototype + BA.1 mRNA",
                             Treatment %in% c("Beta (Pfizer 1)", "Beta mRNA") ~ "Beta mRNA",
                             Treatment %in% c("Beta + Wildtype/Prototype (Pfizer 1)","Beta + Prototype mRNA", "Beta + Wildtype/Prototype (Pfizer 1)") ~ "Prototype + Beta mRNA",
                             Treatment %in% c("Beta (Sanofi)","Beta Protein") ~ "Beta Protein",
                             Treatment %in% c("Beta + Prototype (Sanofi)","Beta + Prototype Protein") ~ "Prototype + Beta Protein",
                             Treatment %in% c("Prototype (Sanofi)","Prototype Protein") ~ "Prototype Protein",
                             Treatment %in% c("Omicron (Pfizer 1)", "Omicron BA.1 mRNA") ~ "Omicron BA.1 mRNA",
                             Treatment == "Beta + Omicron mRNA" ~ "Beta + BA.1 mRNA",
                             Treatment == "Delta + Omicron mRNA" ~ "Delta + BA.1 mRNA",
                             TRUE ~ "Error"),
         Infection = ifelse(`Subject ID` %in% .$`Subject ID`[.$infect_flag == "Y"], "Y", "N"),
         Timepoint = factor(Timepoint, levels = c("1", "15", "71", "90", "180")),
         Booster = factor(Booster, levels = c("Prototype mRNA", "Prototype Protein", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA",
                                              "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein",
                                              "Beta + BA.1 mRNA", "Delta + BA.1 mRNA"))) %>%
  mutate(across(contains(c('Proto+', 'Beta+', 'Delta+', 'BA.1+', 'XBB+')), \(x) replace_na(x,0))) %>%
  select(-contains("Live/IgG/Beta"), -contains("Live/IgG/Proto")) %>%
  group_by(`Subject ID`)%>%
  mutate(primary_vax_series = unique(primary_vax_series)[1],
         infect_baseline = unique(infect_baseline)[1],
         pre_study_booster = unique(pre_study_booster)[1],
         EarlyInfection = case_when(c("1_Y") %in% unique(paste0(Timepoint, "_", infect_flag)) ~ "Y",
                                    c("15_Y") %in% unique(paste0(Timepoint, "_", infect_flag)) ~ "Y",
                                    TRUE ~ "N")) %>%
  ungroup()

flow2$Platform <- str_extract(flow2$Booster, "(mRNA)|(Protein)")
flow2$Immunogen <- case_when(flow2$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                             flow2$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                             flow2$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                             flow2$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                             flow2$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                             flow2$Booster == "Beta + BA.1 mRNA" ~ "Beta + BA.1",
                             flow2$Booster == "Delta + BA.1 mRNA" ~ "Delta + BA.1",
                             TRUE ~ flow2$Booster)

########
#Filter the data
###
#create variables to filter out the data
flow2$infect_flag <- ifelse(is.na(flow2$infect_flag), "0", flow2$infect_flag)
flow2$infect_baseline <- ifelse(is.na(flow2$infect_baseline), "N", flow2$infect_baseline)
flow2$RepeatBarcode <- paste(flow2$`Subject ID`, flow2$Timepoint)
repeats <- flow2$`Subject ID`[duplicated(flow2$RepeatBarcode)]

#filter time!!
flow3 <- flow2 %>%
  filter(!`Subject ID` %in% oosBoost) %>%
  filter(!`Subject ID` %in% c(5053514948, 5352514948)) %>% #these 2 donors were infected d1-d15 but infect_flag was erroneously false
  filter(`Subject ID` != 5752524948, #this donor has extraordinarily few memory B cells, so we should remove them
         !(`Subject ID` == 5349564848 & Timepoint == "15")) #this donor had problems with viability for day 15 sample)
flow3 <- flow3 %>% filter(!(`Subject ID` %in% repeats & Dataset == "Delta Panel"))

#remove infected donors at baseline *for now*
flow3  <- flow3 %>% filter(infect_baseline != "Y")

#remove donors missing data from either day 15 and/or day 0
day0 <- flow3$`Subject ID`[flow3$Timepoint == 1]
day15 <- flow3$`Subject ID`[flow3$Timepoint == 15]
missing <- flow3$`Subject ID`[!(flow3$`Subject ID` %in% day0 & flow3$`Subject ID` %in% day15)]
flow3 <- flow3 %>% filter(!`Subject ID` %in% missing, !is.na(Timepoint), Timepoint != 71) #remove missing donors as well as samples processed with missing timepoints
# `Subject ID` != 5453494948, `Subject ID` != 5151544948) #donors with a *massive* anti-RBD response and a super high baseline

#exclude flow data from delta+ba1 groups
flow3 <- flow3 %>% filter(Booster != "Delta + BA.1 mRNA")
#####

#####
# Creating a cleaner document
intermediateFlow2 <- flow3 %>%
  mutate(`Combined_Proto+/Beta+/BA1+` = rowSums(select(., `Proto+/Beta+/BA1+`,`Proto+/Beta+/BA1+/Delta+`, `Proto+/Beta+/BA1+/XBB+`)),
         `Combined_Proto+/Beta+` = rowSums(select(., `Proto+/Beta+`,`Proto+/Beta+/Delta+`, `Proto+/Beta+/XBB+`)),
         `Combined_Proto+/BA1+` = rowSums(select(., `Proto+/BA1+`,`Proto+/BA1+/Delta+`, `Proto+/BA1+/XBB+`)),
         `Combined_Proto+` = rowSums(select(., `Proto+`,`Proto+/Delta+`, `Proto+/XBB+`)),
         `Combined_Beta+/BA1+` = rowSums(select(., `Beta+/BA+`,`Beta+/BA1+/Delta+`, `Beta+/BA1+/XBB+`)),
         `Combined_Beta+` = rowSums(select(., `Beta+`,`Beta+/Delta+`, `Beta+/XBB+`)),
         `Combined_BA1+` = rowSums(select(., `BA1+`,`BA1+/Delta+`, `BA1+/XBB+`)))

finalFlow2 <- intermediateFlow2 %>%
  select(`Subject ID`, `Specimen ID`, Timepoint, primary_vax_series, infect_baseline, pre_study_booster, Infection, infect_flag, oosboost_flag, Booster, Immunogen, Company, Platform, Stage, Arm,
         Dataset, `Lymphocytes | Count`, `Lymphocytes | Freq. of Total`, `Lymphocytes/Single Cells/CD19 | Count`, `Lymphocytes/Single Cells/CD19 | Freq. of Single Cells`,
         `Live/IgG/IgG++ | Count`, `Live/IgG/IgG++ | Freq. of Live`,
         `Combined_Proto+/Beta+/BA1+`, `Combined_Proto+/Beta+`, `Combined_Proto+/BA1+`, `Combined_Proto+`, `Combined_Beta+/BA1+`, `Combined_Beta+`, `Combined_BA1+`)

write_xlsx(finalFlow2, here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Filtered_COVAILDataset_Infected.xlsx"))