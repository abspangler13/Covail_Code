library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

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
               "Delta + BA.1 mRNA" = "firebrick4",
               "Beta + BA.1 mRNA" = "wheat1")

#let's work on this
immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#068041",
                     "Prototype + Beta" = "#c47002",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f",
                     "Delta + BA.1" = "firebrick4",
                     "Beta + BA.1" = "wheat1")

#load in the flow data run on XBB
metadata <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx"))
metadata2 <- read.csv(here::here("01_raw-data", "FlowData","bcell_unblinding.csv"))

flow <- read_xlsx(here::here("01_raw-data", "FlowData","241004_CombinedXBBAndDeltaDatasets.xlsx")) %>%
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
         Company = str_extract(CompleteBoost, "(Sanofi)|(Pfizer)|(Moderna)")) %>%
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

#add in the sum totals
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
flow$infect_baseline <- ifelse(is.na(flow$infect_baseline), "N", flow$infect_baseline)

#I've noticed that we have a ton of repeat data points where samples were run on Delta then again on XBB- let's remove!
flow$RepeatBarcode <- paste(flow$`Subject ID`, flow$Timepoint)
repeats <- flow$`Subject ID`[duplicated(flow$RepeatBarcode)]

flow <- flow %>% filter(!(`Subject ID` %in% repeats & Dataset == "Delta Panel"))

#write to stats folder
write.csv(flow, here::here("04_Analysis", "data_objects", "paperfigures", "CombinedFlowData_CompleteDeltaData.csv"))

#remove infected donors *for now*
flow  <- flow %>% filter(infect_baseline != "Y", infect_flag != "Y")

#remove donors missing data from either day 15 and/or day 0
day0 <- flow$`Subject ID`[flow$Timepoint == 1]
day15 <- flow$`Subject ID`[flow$Timepoint == 15]
missing <- flow$`Subject ID`[!(flow$`Subject ID` %in% day0 & flow$`Subject ID` %in% day15)]
flow <- flow %>% filter(!`Subject ID` %in% missing, !is.na(Timepoint), Timepoint != 71) #remove missing donors as well as samples processed with missing timepoints
# `Subject ID` != 5453494948, `Subject ID` != 5151544948) #donors with a *massive* anti-RBD response and a super high baseline

#exclude flow data from delta+ba1 groups
flow <- flow %>% filter(Booster != "Delta + BA.1 mRNA")
#####

#####
#Generate numbers for intro graphic
####
stats <- flow %>%
          group_by(Immunogen, Timepoint) %>%
          summarize(n = length(unique(`Subject ID`)))
####
######
#Fig 1b: Total Response boxplot
flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Prototype + Beta", "Beta", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1", "Delta + BA.1"))

#pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1b_AllBoosters_TotalRBDResponse_fix.pdf"), width = 8, height = 3.5, units = "in")
ggplot(flow[flow$TotalRBD != 0,], aes(x = Timepoint, y=TotalRBD))+ #let's not include infection yet
  stat_boxplot(geom= 'errorbar', width = 0.2, lwd = 0.35)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 1, lwd = 0.1)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.25, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4)+
  ylab("Total RBD+ Memory B Cells (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), rows = vars(Platform), axes = "all", labeller = label_wrap_gen(10))+
  #ylim(c(0,10))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1b_AllBoosters_TotalRBDResponse_fix.png"), width = 7, height = 3.1, units = "in", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1b_AllBoosters_TotalRBDResponse_fix.pdf"), width = 7, height = 3.1, units = "in")
dev.off()

#make stats sheet
statistics <-  flow %>% filter(infect_flag == "0" & Timepoint %in% c(1, 15)) %>%
  select(Booster, `Subject ID`, TotalRBD, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD) %>%
  na.omit()

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "TotalRBD_EveryGroup_NAsRemoved.csv"))
#####

#####
#Fig 1c: FC in total RBD
# stats <- flow %>% filter(Infection == "N" & Immunogen %in% c("Prototype", "Prototype + Beta", "Beta")) %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
#   arrange(Timepoint) %>%
#   mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
#   group_by(Booster,Immunogen,Platform, Timepoint) %>%
#   summarize(length = n(),
#             mean = mean(FoldTotalRBD),
#             se = sd(FoldTotalRBD) / sqrt(length)
#   )
# 
# #pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1d_FoldChangeTotalRBDByPlatform.pdf"),width = 3.5, height = 4)
# ggplot(stats, aes(x = Timepoint, y=mean, fill = Booster))+
#   geom_hline(yintercept = 1, linetype = 2, lwd = 0.4)+
#   geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2, lwd=0.35)+
#   geom_line(aes(group = Booster, color = Booster, linetype = Platform))+
#   geom_point(shape = 21, aes(fill = Booster))+
#   ggtitle("Fold Change in Total RBD")+
#   ylab("Fold Change")+
#   xlab("Days Post-Immunization")+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   facet_grid(rows = vars(Immunogen))+
#   theme_classic()+
#   scale_y_continuous(breaks = c(0.8,1,1.2,1.6,2.0), limits = c(0.7,2.25))+
#   guides(fill = "none", color = "none")+
#   theme(
#     plot.title = element_text(size=9, hjust=0.5), 
#     axis.title.y = element_text(size=8),
#     axis.title.x = element_text(size=8),
#     axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
#     axis.text.y = element_text(size=8),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 8),
#     panel.spacing = unit(0.6, "lines"),
#     legend.text = element_text(size = 6),
#     legend.key.size = unit(0.5, 'cm'),
#     legend.title = element_text(size = 7),
#     legend.margin=margin(0,0,0,0),
#     legend.position = "right",
#     axis.line = element_line(size=0.4))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1d_FoldChangeTotalRBDByPlatform.png"),height = 3.8, width = 2.4, units = "in", device = "png", dpi = 600)
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1d_FoldChangeTotalRBDByPlatform.pdf"),height = 3.8, width = 2.4, units = "in", device = "pdf")
# dev.off()

# #Write stats sheet- I'll do this for all groups
# statistics <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
#   arrange(Timepoint) %>%
#   mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
#   select(Booster, `Subject ID`, FoldTotalRBD, Timepoint) %>%
#   pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)
# 
# write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "FoldChangeInTotalRBD_mrnavsprot.csv"))
# #####

#####
#Fig 1d
stats <- flow %>% filter(Infection == "N") %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  filter(FoldTotalRBD < 20) %>%
  group_by(Immunogen, Timepoint) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            sd = sd(FoldTotalRBD)
  )

#pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Omi.png"),width = 3.5, height = 2)
ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1"),], aes(x = Timepoint, y=mean, fill = Immunogen))+ 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen))+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ggtitle("FC in Total RBD")+
  ylab("Fold Change")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  scale_y_continuous(breaks = c(0.3, 1, 2.0, 3), limits = c(0.3, 3))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(plot.title = element_text(size=9, hjust = 0.5), 
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Omi.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Omi.pdf"),width = 3, height = 1.8)
dev.off()

#pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Beta.png"),width = 3.5, height = 2)
ggplot(stats[stats$Immunogen %in% c("Prototype", "Prototype + Beta", "Beta", "Beta + BA.1"),], aes(x = Timepoint, y=mean, fill = Immunogen))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen))+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ggtitle("FC in Total RBD")+
  ylab("Fold Change")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  scale_y_continuous(breaks = c(0.3, 1, 2.0, 3), limits = c(0.3, 3))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(plot.title = element_text(size=9, hjust=0.5), 
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Beta.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Beta.pdf"),width = 3, height = 1.8, units = "in")
dev.off()

#write a file for stats
stats <- flow %>% filter(Infection == "N" & !Platform == "Protein") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  select(Booster, Immunogen, `Subject ID`, Timepoint, FoldTotalRBD) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)

write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "FCTotalRBD_ImmunogenComparison.csv"))
#####