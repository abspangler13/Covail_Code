library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

#####
#set colors
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

#load in the flow data run on XBB
# flow <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx"))
# flow <- flow[!duplicated(flow),] #remove duplicated entries
flow <- read_xlsx(here::here("01_raw-data", "FlowData","241004_CombinedXBBAndDeltaDatasets.xlsx")) %>%
          filter(!is.na(`Specimen ID`)) %>%
          mutate(`Time point Guess` = as.character(`Time point Guess`),
            Timepoint = case_when(timepoint == "Day 1" ~ "1",
                                       timepoint == "Day 15" ~ "15",
                                       timepoint == "Day 57" ~ "57",
                                       timepoint == "Day 71" ~ "71",
                                       timepoint == "Day 91" ~ "90",
                                       timepoint == "Day 181" ~ "180",
                                       TRUE ~ `Time point Guess`)) %>%
          select(-contains("Live/IgG/Beta"), -contains("Live/IgG/Proto"))

#remove out of study boosts/infections
oosBoost <- unique(flow$`Subject ID`[flow$oosboost_flag == "Y" & !is.na(flow$infect_flag)])
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

flow <- flow %>%
        filter(!`Subject ID` %in% oosBoost) %>%
        mutate(
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
        Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
  mutate(across(contains(c('Proto+', 'Beta+', 'Delta+', 'BA.1+', 'XBB+')), \(x) replace_na(x,0))) %>%
  filter(`Subject ID` != 5752524948, #this donor has extraordinarily few memory B cells, so we should remove them
         !(`Subject ID` == 5349564848 & Timepoint == "15")) %>% #this donor had problems with viability for day 15 sample)
  select(-treatment)

#set timepoint as character and then set order by declaring as factor
flow$Timepoint <- factor(flow$Timepoint, levels = c("1", "15", "71", "90", "180"))

#set desired factor order for booster
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype Protein", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA",
                                                "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein",
                                                "Beta + BA.1 mRNA", "Delta + BA.1 mRNA"))

#Sum each of the totals:
#stage 1- we want a comparison between proto and omi-restricted compartments
flow<- flow %>%
              mutate(TotalRBD = rowSums(select(.,contains(c("Proto+", "Beta+", "BA1+")))),
                     ProtoNotBeta = rowSums(select(., contains("Proto+"), -contains("Beta+"))),
                     BetaNotProto = rowSums(select(., contains("Beta+"), -contains("Proto+"))),
                     ProtoBeta = rowSums( select(.,matches("Proto+.+Beta+"))),
                     ProtoNotOmicron = rowSums(select(., contains("Proto+"), -contains("BA+"), -contains("BA1+"))),
                     OmiNotProto = rowSums(select(., contains("BA+"), contains("BA1+"), -contains("Proto+"))),
                     ProtoOmi = rowSums(select(.,matches("Proto+.+BA+"))))

#set faceting variables now since doing so before would mess with column orders
flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            flow$Booster == "Beta + BA.1 mRNA" ~ "Beta + BA.1",
                            flow$Booster == "Delta + BA.1 mRNA" ~ "Delta + BA.1",
                            TRUE ~ flow$Booster)

flow$infect_flag <- ifelse(is.na(flow$infect_flag), "0", flow$infect_flag)

# #####load the old dataset
# old <- read_xlsx(here::here("04_Analysis", "01_raw-data", "FlowData", "Delta panel final dataset.xlsx"))
# 
# #fix existing variables
# old$Timepoint <- factor(as.character(str_remove(old$Timepoint, "Day ")), levels = c("1", "15", "90", "180"))
# old$infect_flag <- as.character(old$infect_flag)
# old$Stage <- as.double(old$Stage)
# old$`Date run` <- as.double(old$`Date run`)
# 
# #let's calculate all of the standard things for this dataset
# old$TotalRBD <- rowSums(old[,c(22:35)]) #I will exclude delta single-positives here as they are not relevant to any of our vaccination groups
# old$ProtoNotBeta <- rowSums(old[,c(28,31, 32, 34)])
# old$BetaNotProto <- rowSums(old[,c(22,24,25,27)])
# old$ProtoBeta <- rowSums(old[,c(29,30,33,35)])
# old$ProtoNotOmicron <- rowSums(old[,c(28,29,32,33)])
# old$OmicronNotPrototype <- rowSums(old[,c(23,24,26,27)])
# old$ProtoOmi <- rowSums(old[,c(30,31,34,35)])
# old$ProtoBetaNotOmi <- rowSums(old[,c(29,33)])
# old$ProtoOmiNotBeta <- rowSums(old[,c(31,34)])
# old$ProtoBetaOmiCrossReactive <- rowSums(old[,c(35,30)])
# old$probeset <- "Old"
# 
# #only keep columns present in flow
# old <- old[,colnames(old) %in% colnames(flow)]
# 
# #merge with flow data
# flow <- bind_rows(flow, old)
# 
# #rerun old variables so that the old data has them
# infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])
# 
# flow <- flow %>%
#   filter(!`Subject ID` %in% oosBoost) %>%
#   mutate(
#     Treatment = ifelse(Treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", Treatment),
#     Booster = case_when(Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
#                         Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA",
#                         Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
#                         Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
#                         Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
#                         Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
#                         Treatment == "Beta (Sanofi)" ~ "Beta Protein",
#                         Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
#                         Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
#                         Treatment == "Omicron (Pfizer 1)" ~ "Omicron BA.1 mRNA",
#                         Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + BA.1 mRNA",
#                         TRUE ~ "Error"
#     ),
#     Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
#   filter(`Subject ID` != 5752524948) #this donor has extraordinarily few memory B cells, so we should remove them
# 
# flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
# flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
#                             flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
#                             flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
#                             flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
#                             flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
#                             TRUE ~ flow$Booster)
# 
# flow$Immunogen <- factor(flow$Immunogen, levels = c( "Omicron BA.1", "Prototype + BA.1", "Prototype", "Prototype + Beta", "Beta"))
# 
# #add the updated data from the unblinding for infected group
# unblinding <- read.csv(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "metadata", "VRC_infections_Details_01MAR2024_Stg1_2_3.csv"))
# 
# flow <- flow %>%
#           mutate(
#             EnrollmentDate = unblinding$ENRLDATE[match(`Subject ID`, unblinding$Subject_ID)],
#             InfectionDate = unblinding$datinf[match(`Subject ID`, unblinding$Subject_ID)],
#             InfectionLineage = unblinding$PANGOLIN_SWAB1[match(`Subject ID`, unblinding$Subject_ID)],
#             InfectionTrunc = unblinding$TRUNCLIN_SWAB1[match(`Subject ID`, unblinding$Subject_ID)],
#             InfectionIndentMethod = unblinding$source[match(`Subject ID`, unblinding$Subject_ID)]
#           ) %>%
#         filter(`Subject ID` != 5450574848) #remove a specific donor for whom we don't have day 1 data for (but do have days 15 and 180)
# #the donor removed above was not included in statistics

#I've noticed that we have a ton of repeat data points where samples were run on Delta then again on XBB- let's remove!
flow$RepeatBarcode <- paste(flow$`Subject ID`, flow$Timepoint)
repeats <- flow$`Subject ID`[duplicated(flow$RepeatBarcode)]

flow <- flow %>% filter(!(`Subject ID` %in% repeats & Dataset == "Delta Panel"))

#write to stats folder
write.csv(flow, here::here("04_Analysis", "data_objects", "paperfigures", "CombinedFlowData_CompleteDeltaData.csv"))

#remove infected donors *for now*
flow  <- flow %>% filter(!Booster %in% c("Beta + BA.1 mRNA", "Delta + BA.1 mRNA"),
                          Infection != "Y",
                          infect_baseline != "Y" & !is.na(infect_baseline))
  
#remove donors missing data from either day 15 and/or day 0
day0 <- flow$`Subject ID`[flow$Timepoint == 1]
missing <- flow$`Subject ID`[flow$`Subject ID`

#####

#####
#Generate numbers
####
stats <- flow %>%
          group_by(Immunogen, Timepoint) %>%
          summarize(n = length(unique(`Subject ID`)))

####

######
#Fig 1b: Total Response- do we see a specific group dominating in total responses?
flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Prototype + Beta", "Beta", "Prototype + BA.1", "Omicron BA.1"))

ggplot(flow[flow$TotalRBD != 0,], aes(x = Timepoint, y=TotalRBD))+ #let's not include infection yet
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1b_AllBoosters_TotalRBDResponse_fix.png"),width = 7.2, height = 3.5, units = "in", device = "png", dpi = 600)
dev.off()

#make stats sheet
statistics <-  flow %>% filter(infect_flag == "0" & Timepoint %in% c(1, 15)) %>%
  select(Booster, `Subject ID`, TotalRBD, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD) %>%
  na.omit()

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 1", "TotalRBD_EveryGroup_NAsRemoved.csv"))
#####

#####
#Fig 1c: Total Response- do we see a specific group dominating in total responses by fold change?
stats <- flow %>% filter(Infection == "N" & !Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster,Immunogen,Platform, Timepoint) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

ggplot(stats, aes(x = Timepoint, y=mean, fill = Booster))+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1d_FoldChangeTotalRBDByPlatform.png"),width = 3.5, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#Write stats sheet- I'll do this for all groups
statistics <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  select(Booster, `Subject ID`, FoldTotalRBD, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 1", "FoldChangeInTotalRBD_AllGroups.csv"))
#####

#####
#Fig 1d
stats <- flow %>% filter(Infection == "N" & !Platform == "Protein") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster,Immunogen,Platform, Timepoint) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + BA.1", "Omicron BA.1"),], aes(x = Timepoint, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Omi.png"),width = 3.5, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + Beta", "Beta"),], aes(x = Timepoint, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Beta.png"),width = 3.5, height = 2, units = "in", device = "png", dpi = 600)
dev.off()
#####