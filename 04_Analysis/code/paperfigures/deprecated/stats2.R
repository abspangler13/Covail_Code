library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(rstatix)
library(gridExtra)
library(readxl)
library(writexl)
library(openxlsx)

#load in the flow data
flow <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx"))
flow <- flow[!duplicated(flow),] #remove duplicated entries

#add in booster and remove anyone who has had an out-of-study boost
oosBoost <- unique(flow$`Subject ID`[flow$oosboost_flag == "Y"])
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

flow <- flow %>%
  filter(!`Subject ID` %in% oosBoost) %>%
  mutate(
    Treatment = ifelse(Treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", Treatment),
    Booster = case_when(Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                        Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA",
                        Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                        Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                        Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                        Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                        Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                        Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                        Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                        Treatment == "Omicron (Pfizer 1)" ~ "Omicron BA.1 mRNA",
                        Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + BA.1 mRNA",
                        TRUE ~ "Error"
    ),
    Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
  filter(`Subject ID` != 5752524948, #this donor has extraordinarily few memory B cells, so we should remove them
         !(`Subject ID` == 5349564848 & `Time point Guess` == "15")) #this donor had problems with viability for day 15 sample)

#set timepoint as character and then set order by declaring as factor
flow$`Time point Guess` <- as.character(flow$`Time point Guess`)
flow$`Time point Guess` <- factor(flow$`Time point Guess`, levels = c("1", "15", "90", "180"))

#set desired factor order for booster
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype Protein", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA",
                                                "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein"))

#some funky stuff happened when Flavio exported from FlowJo, so we'll need to divide all of the values by 100 from the Pfizer groups
for(i in 1:nrow(flow)){
  if(flow$Treatment[i] %in% c("Omicron (Pfizer 1)", "Omicron + Wildtype/Prototype (Pfizer 1)")){
    for(j in c(seq(from = 24, to = 54, by =2))){
      flow[i,j] <- flow[i,j] / 100
    }
  }
} #wow i slayed that

#Prepare probe metrics
#Proto+Beta+
flow$`Proto+/Beta+/BA1+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+/BA1+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1 | Freq. of IgG` * 100
flow$`Proto+/Beta+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta | Freq. of IgG` * 100

#Proto Only
flow$`Proto+` <- flow$`Live/IgG/Proto+/Proto | Freq. of IgG` * 100
flow$`Proto+/BA1+/XBB+` <- flow$`Live/IgG/Proto+/Proto-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/BA1+` <- flow$`Live/IgG/Proto+/Proto-BA1 | Freq. of IgG` * 100
flow$`Proto+/XBB+` <- flow$`Live/IgG/Proto+/Proto-XBB | Freq. of IgG` * 100

#Beta Only
flow$`Beta+` <- flow$`Live/IgG/Beta+/Beta | Freq. of IgG` * 100
flow$`Beta+/BA1+/XBB+` <- flow$`Live/IgG/Beta+/Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Beta+/BA+` <- flow$`Live/IgG/Beta+/Beta-BA1 | Freq. of IgG` * 100
flow$`Beta+/XBB+` <- flow$`Live/IgG/Beta+/Beta-XBB | Freq. of IgG` * 100

#the m'crons
flow$`BA1+` <- flow$`Live/IgG/Proto neg Beta neg/BA1 | Freq. of IgG` * 100
flow$`XBB+` <- flow$`Live/IgG/Proto neg Beta neg/XBB | Freq. of IgG` * 100
flow$`BA1+/XBB+` <- flow$`Live/IgG/Proto neg Beta neg/BA1-XBB | Freq. of IgG` * 100

#Sum each of the totals:
#stage 1- we want a comparison between proto and omi-restricted compartments
flow$TotalRBD <- rowSums(flow[,c(57:69, 71)])
flow$ProtoNotBeta <- rowSums(flow[,c(61:64)])
flow$BetaNotProto <- rowSums(flow[,c(65:68)])
flow$ProtoBeta <- rowSums(flow[,c(57:60)])
flow$ProtoNotOmicron <- rowSums(flow[,c(59, 60, 61, 64)])
flow$OmicronNotPrototype <- rowSums(flow[,c(66, 67, 69, 71)])
flow$ProtoOmi <- rowSums(flow[,c(57, 58, 62, 63)])
flow$ProtoBetaNotOmi <- rowSums(flow[,c(59, 60)])
flow$ProtoOmiNotBeta <- rowSums(flow[,c(62, 63)])
flow$ProtoBetaOmiCrossReactive <- rowSums(flow[,c(57, 58)])
flow$ProtoBetaBA1NotXBB <- rowSums(flow[,c(58)])
flow$ProtoBetaBA1XBB <- rowSums(flow[,c(57)])
flow$ProtoBA1NotXBB <- rowSums(flow[,c(58, 63)])
flow$ProtoBA1XBB <- rowSums(flow[,c(57, 62)])

#set faceting variables now since doing so before would mess with column orders
flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster)

flow$infect_flag <- ifelse(is.na(flow$infect_flag), "0", flow$infect_flag)

#we will be merging our data with an earlier dataset run with a different probe panel (Beta on diff fluor, Delta variant instead of XBB)
#let's create a variable to differentiate the two
flow$probeset <- "New"

#####load the old dataset
old <- read_xlsx(here::here("01_raw-data", "FlowData","Stage23DeltaPaneltoincludeD0D15.xlsx"), sheet = "newdata")
old <- old[!duplicated(old),] #remove duplicated entries

#fix existing variables
old$`Time point Guess` <- factor(as.character(str_remove(old$`Time point Guess`, "Day ")), levels = c("1", "15", "90", "180"))
old$infect_flag <- as.character(old$infect_flag)
old$Stage <- as.double(old$Stage)
old$`Date run` <- as.double(old$`Date run`)

#let's calculate all of the standard things for this dataset
old$TotalRBD <- rowSums(old[,c(22:35)]) #I will exclude delta single-positives here as they are not relevant to any of our vaccination groups
old$ProtoNotBeta <- rowSums(old[,c(28,31, 32, 34)])
old$BetaNotProto <- rowSums(old[,c(22,24,25,27)])
old$ProtoBeta <- rowSums(old[,c(29,30,33,35)])
old$ProtoNotOmicron <- rowSums(old[,c(28,29,32,33)])
old$OmicronNotPrototype <- rowSums(old[,c(23,24,26,27)])
old$ProtoOmi <- rowSums(old[,c(30,31,34,35)])
old$ProtoBetaNotOmi <- rowSums(old[,c(29,33)])
old$ProtoOmiNotBeta <- rowSums(old[,c(31,34)])
old$ProtoBetaOmiCrossReactive <- rowSums(old[,c(35,30)])
old$probeset <- "Old"

#only keep columns present in flow
old <- old[,colnames(old) %in% colnames(flow)]

#merge with flow data
flow <- bind_rows(flow, old)

#rerun old variables so that the old data has them
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

flow <- flow %>%
  filter(!`Subject ID` %in% oosBoost) %>%
  mutate(
    Treatment = ifelse(Treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", Treatment),
    Booster = case_when(Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                        Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA",
                        Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                        Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                        Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                        Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                        Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                        Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                        Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                        Treatment == "Omicron (Pfizer 1)" ~ "Omicron BA.1 mRNA",
                        Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + BA.1 mRNA",
                        TRUE ~ "Error"
    ),
    Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))

flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster)

flow$Immunogen <- factor(flow$Immunogen, levels = c( "Omicron BA.1", "Prototype + BA.1", "Prototype", "Prototype + Beta", "Beta"))
#####

#####Figure 2a: Cross reactives over time
stats <- flow %>% filter(infect_flag == "0" & Booster %in% c("Prototype mRNA", "Prototype Protein", "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein")) %>%
  select(Treatment, `Subject ID`,probeset, ProtoBeta, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = ProtoBeta) %>%
  mutate(AnyNAs = rowSums(is.na(.)) >= 1)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2a_ProtoBetaCrossReactive_OverTime.xlsx"))

stats <- flow %>% filter(infect_flag == "0" & Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, ProtoOmi, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = ProtoOmi) %>%
  mutate(AnyNAs = rowSums(is.na(.)) >= 1)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2a_ProtoOmiCrossReactive_OverTime.xlsx"))
#####

#####Figure 2b: Prototype specific over time
stats <- flow %>% filter(infect_flag == "0" & Booster %in% c("Prototype mRNA", "Prototype Protein", "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein")) %>%
  select(Booster, `Subject ID`,probeset, ProtoNotBeta, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = ProtoNotBeta) %>%
  mutate(AnyNAs = rowSums(is.na(.)) >= 1)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2b_ProtoNotBeta_OverTime.xlsx"))

stats <- flow %>% filter(infect_flag == "0" & Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>%
  select(Booster, `Subject ID`,probeset, ProtoNotOmicron, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = ProtoNotOmicron) %>%
  mutate(AnyNAs = rowSums(is.na(.)) >= 1)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2b_ProtoNotOmicronBA1_OverTime.xlsx"))
#####

#####Figure 2d: ratio of cross:prototype-specific
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

stats <- flow %>% filter(infect_flag == "0" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Prototype Protein", "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein")) %>%
  select(Treatment, `Subject ID`,probeset, RatioCrossBeta, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = RatioCrossBeta)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2d_RatioProtoBetaToProtoNotBeta_OverTime.xlsx"))

stats <- flow %>% filter(infect_flag == "0" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, RatioCrossOmi, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = RatioCrossOmi)

write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure2", "Figure2d_RatioProtoOmiToProtoNotOmi_OverTime.xlsx"))

# do fold change
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA")) %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  RatioCrossOmi / RatioCrossOmi [1]) %>%
  select(Treatment, `Subject ID`,probeset, Fold, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = Fold)
write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure1", "Figure2f_FOLDCHANGE_RatioProtoOmiToProtoNotOmi_OverTime.xlsx"))


stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein")) %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  RatioCrossBeta / RatioCrossBeta [1]) %>%
  select(Treatment, `Subject ID`,probeset, Fold, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = Fold)
write_xlsx(stats, here::here("04_Analysis","plots","paperfigures", "stats", "Figure1", "Figure2f_FOLDCHANGE_RatioProtoBetaToProtoNotBeta_OverTime.xlsx"))
#####