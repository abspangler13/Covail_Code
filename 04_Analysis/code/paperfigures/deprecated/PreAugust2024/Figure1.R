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
library(writexl)

#####
#set colors
#we'll be coloring by immunogen rather than platform now
# allColors <- c("Omicron BA.1 mRNA" = "#FF1F5B", 
#                  "Prototype + BA.1 mRNA" = "#009ADE",
#                  "Prototype mRNA" = "#FFC61E",
#                   "Prototype Protein" = "#ff6e54",
#                   "Prototype + Beta mRNA" = "#dd5182",
#                   "Prototype + Beta Protein" = "#955196",
#                   "Beta mRNA" = "#6573bf",
#                   "Beta Protein" = "#005e8a")

#based on Zissou1 color palette from wes anderson
allColors <- c("Prototype mRNA" = "#045275",
               "Prototype Protein" = "#0a86bf",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#f78c00",
               "Prototype + BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f")

immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#7CCBA2",
                     "Prototype + Beta" = "#FCD39C",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f")

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
    Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
  filter(`Subject ID` != 5752524948) #this donor has extraordinarily few memory B cells, so we should remove them

flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster)

flow$Immunogen <- factor(flow$Immunogen, levels = c( "Omicron BA.1", "Prototype + BA.1", "Prototype", "Prototype + Beta", "Beta"))

#add the updated data from the unblinding for infected group
unblinding <- read.csv(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "metadata", "VRC_infections_Details_01MAR2024_Stg1_2_3.csv"))

flow <- flow %>%
          mutate(
            EnrollmentDate = unblinding$ENRLDATE[match(`Subject ID`, unblinding$Subject_ID)],
            InfectionDate = unblinding$datinf[match(`Subject ID`, unblinding$Subject_ID)],
            InfectionLineage = unblinding$PANGOLIN_SWAB1[match(`Subject ID`, unblinding$Subject_ID)],
            InfectionTrunc = unblinding$TRUNCLIN_SWAB1[match(`Subject ID`, unblinding$Subject_ID)],
            InfectionIndentMethod = unblinding$source[match(`Subject ID`, unblinding$Subject_ID)]
          )

#write to stats folder
writexl::write_xlsx(flow, here::here("04_Analysis", "data_objects", "paperfigures", "CombinedFlowData.xlsx"))
#####

######
#Fig 1b: Total Response- do we see a specific group dominating in total responses?
ggplot(flow[flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=TotalRBD))+ #let's not include infection yet
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd = 0.4)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  ylab("Total RBD+ Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), rows = vars(Platform), axes = "all", labeller = label_wrap_gen(15))+
  #ylim(c(0,10))+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1b_AllBoosters_TotalRBDResponse.png"),width = 7.4, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#make stats sheet
statistics <-  flow %>% filter(infect_flag == "0" & `Time point Guess` %in% c(1, 15)) %>%
  select(Booster, `Subject ID`, TotalRBD, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = TotalRBD) %>%
  na.omit()

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 1", "TotalRBD_EveryGroup_NAsRemoved.csv"))
#####

#####
#Fig 1c: Total Response- do we see a specific group dominating in total responses by fold change?
stats <- flow %>% filter(Infection == "N" & !Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster,Immunogen,Platform, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster, linetype = Platform))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Total RBD")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(rows = vars(Immunogen))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=9), 
    axis.title.y = element_text(size=8),
    axis.title.x = element_text(size=8),
    axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.6, "lines"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.1, 'cm'),
    legend.title = element_text(size = 7),
    legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1d_FoldChangeTotalRBDByPlatform.png"),width = 3.5, height = 4.5, units = "in", device = "png", dpi = 600)
dev.off()

#Write stats sheet- I'll do this for all groups
statistics <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  select(Booster, `Subject ID`, FoldTotalRBD, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = FoldTotalRBD)

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 1", "FoldChangeInTotalRBD_AllGroups.csv"))
#####

#####
#Fig 1d
stats <- flow %>% filter(Infection == "N" & !Platform == "Protein") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster,Immunogen,Platform, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + BA.1", "Omicron BA.1"),], aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Total RBD")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(plot.title = element_text(size=9), 
    axis.title.y = element_text(size=8),
    axis.title.x = element_text(size=8),
    axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8),
    panel.spacing = unit(0.6, "lines"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.1, 'cm'),
    legend.title = element_text(size = 7),
    legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Omi.png"),width = 3.5, height = 2.25, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + Beta", "Beta"),], aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Total RBD")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(plot.title = element_text(size=9), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.6, "lines"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 7),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Beta.png"),width = 3.5, height = 2.25, units = "in", device = "png", dpi = 600)
dev.off()
#####