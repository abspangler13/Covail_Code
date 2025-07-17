library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(rstatix)
library(gridExtra)
library(thematic)
library(readxl)
library(viridisLite)
library(stringdist)
library(corrplot)
library(ggpubr)
library(ggprism)
library(wesanderson)
library(rstatix)
library(gridExtra)
library(openxlsx)

#####
#set colors
allColors <- c("Prototype mRNA" = "#045275",
               "Prototype Protein" = "#0a86bf",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#f78c00",
               "Prototype + BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f",
               "Beta + BA.1 mRNA" = "firebrick",
               "Delta + BA.1 mRNA" = "goldenrod1")

immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#7CCBA2",
                     "Prototype + Beta" = "#FCD39C",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f")

#####load the old dataset
old <- read_xlsx(here::here("01_raw-data", "FlowData","Stage23DeltaPaneltoincludeD0D15.xlsx"), sheet = "newdata")
old <- old[!duplicated(old),] #remove duplicated entries

#fix existing variables
old$Timepoint <- old$`Time point Guess`
old$infect_flag <- as.character(old$infect_flag)
old$Stage <- as.double(old$Stage)
old$`Date run` <- as.double(old$`Date run`)
old$Subject <- old$`Subject ID`

#let's calculate all of the standard things for this dataset
old$probeset <- "Old"

#####
#Add in really old flow data using delta as a probe with all the weird vaccination groups
reallyOld <- read_xlsx(here::here("01_raw-data", "FlowData", "Stage1_Analysis_SFA_olddatamaybe.xlsx"), sheet = "Sheet5")
reallyOldMetadata <- read.csv(here::here("01_raw-data", "FlowData", "bcell_unblinding.csv"))
reallyOldMetadata2 <- read_xlsx(here::here("01_raw-data", "FlowData", "Stage1_Analysis_SFA_olddatamaybe.xlsx"), sheet = "count")

reallyOld <- reallyOld %>% mutate(Subject = reallyOldMetadata2$`Subject#`[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
                                  Timepoint = reallyOldMetadata$timepoint[match(`Sample:`, reallyOldMetadata$sn)],
                                  TimepointCheck = reallyOldMetadata2$Timepoint[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
                                  Treatment = reallyOldMetadata$treatment[match(`Sample:`, reallyOldMetadata$sn)],
                                  TreatmentCheck = reallyOldMetadata2$Group[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
                                  infect_baseline = reallyOldMetadata$infect_baseline[match(`Sample:`, reallyOldMetadata$sn)],
                                  probeset = "Really Old") %>%
  filter(!is.na(Subject)) %>%
  mutate(Treatment = case_when(TreatmentCheck == "T2" ~ "1 Dose Omicron + Prototype (Moderna)",
                               TreatmentCheck == "T3" ~ "1 Dose Omicron (Moderna)",
                               TreatmentCheck == "T4" ~ "1 Dose Beta + Omicron (Moderna)",
                               TreatmentCheck == "T5" ~ "1 Dose Prototype (Moderna)",
                               TRUE ~ Treatment)) %>%
  rename(`Beta+` = `Live/IgG/IgG++/Omi-Proto neg/Beta | Freq. of IgG++`,
         `Delta+` = `Live/IgG/IgG++/Omi-Proto neg/Delta | Freq. of IgG++`,
         `BA1+`=`Live/IgG/IgG++/Omi+ Prot-/Omi | Freq. of IgG++`,
         `Proto+`=`Live/IgG/IgG++/Proto+ Omi-/Proto | Freq. of IgG++`,
         `Proto+/Beta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Beta | Freq. of IgG++`,
         `Proto+/Delta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Delta | Freq. of IgG++`,
         `Proto+/Beta+/Delta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Delta-Beta | Freq. of IgG++`,
         `Proto+/BA1+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi | Freq. of IgG++`,
         `Proto+/Beta+/BA1+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Beta | Freq. of IgG++`,
         `Proto+/Beta+/BA1+/Delta+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Beta-Delta | Freq. of IgG++`,
         `Proto+/BA1+/Delta+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Delta | Freq. of IgG++`,
         `Beta+/Delta+` = `Live/IgG/IgG++/Omi-Proto neg/Delta-Beta | Freq. of IgG++`,
         `Beta+/BA1+/Delta+`=`Live/IgG/IgG++/Omi+ Prot-/Omi-Beta-Delta | Freq. of IgG++`,
         `Beta+/BA1+`=`Live/IgG/IgG++/Omi+ Prot-/Omi-Beta | Freq. of IgG++`,
         `BA1+/Delta+` = `Live/IgG/IgG++/Omi+ Prot-/Omi-Delta | Freq. of IgG++`,
         SN = `Sample:`)

#multiply by 100 to do percentage of IgG
# check <- reallyOld %>% mutate(across(`Live/IgG/IgG++/Omi-Proto neg/Beta | Freq. of IgG++`:`Live/IgG/IgG++/Omi+ Prot-/Omi-Delta | Freq. of IgG++`, ~ .* 100))
#commenting out above because flow data looks like it's already been converted to percentage of IGG
df <- reallyOld[,colnames(reallyOld) %in% colnames(old)]
df <- rbind(df, old[,colnames(old) %in% colnames(df)])

#there are a bunch of missing values for the date and the timepoint- let's get them fixed
# df <- df %>%
#       group_by(Subject) %>%
#       mutate(Timepoint = case_when(n() == 1 ~ "Only One Sample Run",
#                                    "Day 1" %in% unique(Timepoint) & is.na(Timepoint) ~ "Day 15",
#                                    "Day 15" %in% unique(Timepoint) & is.na(Timepoint) ~ "Day 1",
#                                    TRUE ~ Timepoint
#                                    ))
#commented out the above for now. There are some subjects who had 1 timepoint run repeatedly. I have no idea
#what to do with these because they could be mislabelled and run at day 57. I'm going to just exclude them for now

#exclude the subjects with 1 or 3 runs
df <- df %>%
  group_by(Subject) %>%
  filter(n() != 3 & n() != 1)

#now we can calculate different metrics
df <- df %>%
  mutate(TotalRBD = rowSums(across(`Beta+`:`BA1+/Delta+`)),
         VariantNotProto = rowSums(select(across(`Beta+`:`BA1+/Delta+`), !starts_with("Proto+"))),
         NonProtoSpecific = rowSums(select(across(`Beta+`:`BA1+/Delta+`), -c(`Proto+`))),
         OfficialBooster = case_when(Treatment == "1 Dose Beta + Omicron (Moderna)" ~ "Beta + BA.1 mRNA",
                                     Treatment == "1 Dose Delta + Omicron (Moderna)" ~ "Delta + BA.1 mRNA",
                                     Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                                     Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                                     Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                                     Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                                     Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                                     Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                                     Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                                     Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                                     TRUE ~ "We have a problem"),
         Immunogen = str_extract(OfficialBooster, ".+(?=( mRNA)|( Protein))"),
         Platform = str_extract(OfficialBooster, "(mRNA)|(Protein)"))

#newer data
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

#fix the variables
flow$SN <- as.character(flow$SN)

#bind the datasets
mergedData <- bind_rows(flow, df) %>%
    mutate(CombinedInfectionFlag = case_when(!(is.na(Infection)) & Infection == "Y" ~ "Y",
                                 !(is.na(infect_baseline)) & infect_flag == "Y" ~ "Y",
                                 !(is.na(infect_baseline)) & infect_baseline == "Y" ~ "Y",
                                 TRUE ~ "N"),
           CombinedBooster = case_when(!is.na(OfficialBooster) ~ OfficialBooster,
                                       !is.na(Booster) ~ Booster,
                                       TRUE ~ "Check again"),
           CombinedTime = case_when(!is.na(Timepoint) ~ str_remove(Timepoint, "Day "),
                                    TRUE ~ `Time point Guess`),
           CombinedSubject = case_when(!is.na(Subject) ~ Subject,
                                       TRUE ~ `Subject ID`)) %>%
    filter(!duplicated(SN), CombinedInfectionFlag == "N")
#####

##
#plot the differences in total RBD by dataset overall.
ggplot(mergedData, aes(x = CombinedTime, y=TotalRBD))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = CombinedSubject, color = CombinedBooster), alpha = 0.3, lwd = 0.4)+
  geom_boxplot(aes(fill = CombinedBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  ylab("Total RBD+ Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(rows = vars(CombinedBooster), cols = vars(probeset), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "totalRBD.png"),width = 4, height = 12, units = "in", device = "png", dpi = 600)
dev.off()


#side by side
ggplot(mergedData, aes(x = CombinedTime, y=TotalRBD))+
  geom_boxplot(aes(fill = probeset), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  ylab("Total RBD+ Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(cols = vars(CombinedBooster), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "totalRBD_sidebyside.png"),width = 12, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#prototype only
ggplot(mergedData, aes(x = CombinedTime, y=`Proto+`))+
  geom_boxplot(aes(fill = probeset), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  ylab("Prototype+ only Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(cols = vars(CombinedBooster), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  #scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "prototypesinglepositives_sidebyside.png"),width = 12, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#There are a lot more prototype single-positives in the newer panel, which makes sense because there will be more
#proto/delta+ B cells than proto+/XBB+, so proto+ should be more common in the new panel. However, if we add a context (i.e. proto/beta), does it differ?
mergedData <- mergedData %>% mutate(
  ProtoNotBeta = rowSums(select(.,`Proto+`, `Proto+/BA1+/XBB+`, `Proto+/BA1+`, `Proto+/XBB+`, `Proto+/Delta+`, `Proto+/BA1+/Delta+`), na.rm=TRUE),
  ProtoNotBA1 = rowSums(select(., `Proto+`, `Proto+/Beta+/XBB+`, `Proto+/Beta+`, `Proto+/XBB+`, `Proto+/Delta+`, `Proto+/Beta+/Delta+`), na.rm=TRUE),
  SubjectIDCombined = case_when(probeset == "New" ~ `Subject ID`,
                                TRUE ~ Subject),
  SubjectTime = paste(mergedData$SubjectIDCombined, mergedData$CombinedTime))

#interestingly, it seems we have some stragglers from the really old dataset that made it into the new dataset (i.e. we reran them with the XBB panel)
#for now, we'll exclude them because our plots will look ugly otherwise, but it's a good internal control
filteredMergedData <- mergedData %>% filter(!(duplicated(SubjectTime) | duplicated(SubjectTime, fromLast = TRUE)))

ggplot(filteredMergedData, aes(x = CombinedTime, y=ProtoNotBeta))+
  geom_boxplot(aes(fill = probeset), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  #geom_point(shape =21, aes(fill = probeset))+
  #geom_line(aes(color = probeset, group = SubjectIDCombined))+
  ylab("Prototype+/Beta- Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(cols = vars(CombinedBooster), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "prototypepositive_betanegative_sidebyside.png"),width = 12, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(filteredMergedData, aes(x = CombinedTime, y=ProtoNotBA1))+
  geom_boxplot(aes(fill = probeset), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  ylab("Prototype+/BA1- Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(cols = vars(CombinedBooster), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "prototypepositive_ba1negative_sidebyside.png"),width = 12, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#let's plot the people shared between the really old dataset and the new one
repeats <- mergedData %>% filter((duplicated(SubjectTime) | duplicated(SubjectTime, fromLast = TRUE)))

ggplot(repeats, aes(x = CombinedTime, y=ProtoNotBA1))+
  #geom_boxplot(aes(fill = probeset), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_point(shape =21, aes(fill = as.character(SubjectIDCombined)))+
  geom_line(aes(color = as.character(SubjectIDCombined), group = SubjectIDCombined))+
  ylab("Prototype+/BA.1- Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  facet_grid(cols = vars(CombinedBooster), rows = vars(probeset), axes = "all", labeller = label_wrap_gen(10))+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  #scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "flowqc", "prototypepositive_ba.1negative_controls.png"),width = 5, height = 4, units = "in", device = "png", dpi = 600)
dev.off()
