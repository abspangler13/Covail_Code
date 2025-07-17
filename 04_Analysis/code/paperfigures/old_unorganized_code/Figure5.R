library(ggplot2)
library(dplyr)
library(here)
library(Seurat)
library(readxl)
library(tidyseurat)
library(stringr)
library(alakazam)
library(writexl)
library(gridExtra)
library(ggalluvial)

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

#table of sample numbers
stats <- flowRaw %>%
  group_by(Immunogen, Platform, Timepoint) %>%
  summarize(n = length(unique(`Subject ID`)))

#####

#####
#Do we see changes in Total RBD over time?
infectedFlow <- flow  %>%
      mutate(TimeInf = paste0(Timepoint, "_", infect_flag)) %>%
      group_by(`Subject ID`) %>%
      mutate(Range = case_when("1_Y" %in% TimeInf ~ "Uh oh",
                               "15_Y" %in% TimeInf ~ "Uh oh",
                               "90_Y" %in% TimeInf ~ "Infected Days 15-90",
                               "180_Y" %in% TimeInf ~ "Infected Days 90-180",
                               Infection == "N" ~ "Uninfected")) %>%
      mutate(Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180"))) %>%
      filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
      filter((`Subject ID` %in% seuObj$Subject) | (Range %in% c("Infected Days 15-90", "Infected Days 90-180")))

#by flow do we see changes in total RBD over time between the two groups?
infectedFlow %>% mutate(Booster = factor(Booster, levels = c("Prototype mRNA", "Prototype/BA.1 mRNA", "Omicron BA.1 mRNA"))) %>%
ggplot(aes(x = Timepoint, y = TotalRBD))+
  geom_line(aes(color = Range, group = `Subject ID`), alpha= 0.5)+
  #geom_boxplot(aes(fill = Booster), outlier.shape = 21)+
  ylab("Total RBD+ (% of IgG)")+
  geom_point(aes(fill = Range), shape =21, size = 1)+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  facet_grid(cols = vars(Range), axes = "all")+
  theme_classic()+
  theme(text = element_text(size = 8),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "Flow_TotalRBD_Infected.png"), width = 3.7, height = 1.9, dpi = 1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "Flow_TotalRBD_Infected.svg"), width = 3.7, height = 1.9)

#write stats sheet
infectedFlow %>%
  select(`Subject ID`, Timepoint, Range, TotalRBD) %>%
  group_by(`Subject ID`, Range) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD) %>% select(`Subject ID`, Range, `1`, `15`, `90`, `180`) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "Flow_TotalRBD_Infected.xlsx"))
#####

#####
# #look at fold change
# allFlow <- flow %>%
#   mutate(TimeInf = paste0(Timepoint, "_", infect_flag)) %>%
#   group_by(`Subject ID`) %>%
#   mutate(Range = case_when("1_Y" %in% TimeInf ~ "Uh oh",
#                            "15_Y" %in% TimeInf ~ "Uh oh",
#                            "90_Y" %in% TimeInf ~ "Infected Days 15-90",
#                            "180_Y" %in% TimeInf ~ "Infected Days 90-180",
#                            TRUE ~ "Uninfected"),
#          Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180"))) %>%
#   filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
#   mutate(Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180"))) %>%
#   filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
#   filter((`Subject ID` %in% seuObj$Subject) | (Range %in% c("Infected Days 15-90", "Infected Days")))
# 
# calcs <- allFlow %>%
#             group_by(Range, `Subject ID`) %>%
#             arrange(Timepoint)%>%
#             mutate(Fold = ProtoOmi / ProtoOmi[1]) %>%
#             group_by(Range, Timepoint) %>%
#             summarize(n = n(),
#                       mean = mean(Fold),
#                       sd = sd(Fold)) %>%
#             mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("Fold Change")+
#   ggtitle("Prototype+BA.1+")+
#   ylim(0.5, 7)+
#   geom_hline(yintercept = 1, linetype = 2, linewidth= 0.8)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.position = "none",
#         plot.title = element_text(hjust = 0.5, vjust = 0.5))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow.png"), width = 1.3, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow.svg"), width = 1.3, height =1.8)
# 
# 
# #save sheet
# allFlow %>%
#   group_by(Range, `Subject ID`) %>%
#   arrange(Timepoint)%>%
#   mutate(Fold = ProtoOmi / ProtoOmi[1]) %>%
#   select(Range, `Subject ID`, Timepoint, Fold) %>% group_by(Timepoint, `Subject ID`) %>%
#   pivot_wider(names_from = Range, values_from = Fold) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "CrossReactive_FoldChange_Flow.xlsx"))
# 
# 
# ###################omicron specific?
# calcs <- allFlow %>%
#   group_by(Range, `Subject ID`) %>%
#   arrange(Timepoint)%>%
#   mutate(Fold = OmiNotProto / OmiNotProto[1]) %>%
#   group_by(Range, Timepoint) %>%
#   summarize(n = n(),
#             mean = mean(Fold),
#             sd = sd(Fold)) %>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("Fold Change")+
#   ggtitle("Prototype-BA.1+")+
#   ylim(0.5, 7)+
#   geom_hline(yintercept = 1, linetype = 2, linewidth= 0.8)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0.5),
#         legend.key.size = unit(0.4, "lines"))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow.png"), width = 2.5, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow.svg"), width = 2.5, height =1.8)
# 
# allFlow %>% 
#   group_by(Range, `Subject ID`) %>%
#   arrange(Timepoint)%>%
#   mutate(Fold = OmiNotProto / OmiNotProto[1]) %>%
#   select(Range, `Subject ID`, Timepoint, Fold) %>% group_by(Timepoint, `Subject ID`) %>%
#   pivot_wider(names_from = Range, values_from = Fold) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "BA1Only_FoldChange_Flow.xlsx"))
# 
# ###################proto-specific?
# calcs <- allFlow %>%
#   group_by(Range, `Subject ID`) %>%
#   arrange(Timepoint)%>%
#   mutate(Fold = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
#   group_by(Range, Timepoint) %>%
#   summarize(n = n(),
#             mean = mean(Fold),
#             sd = sd(Fold)) %>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("Fold Change")+
#   ggtitle("Prototype+BA.1-")+
#   ylim(0.5, 7)+
#   geom_hline(yintercept = 1, linetype = 2, linewidth= 0.8)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0.5),
#         legend.key.size = unit(0.4, "lines"))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow.png"), width = 2.5, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow.svg"), width = 2.5, height =1.8)
# 
# allFlow %>% 
#   group_by(Range, `Subject ID`) %>%
#   arrange(Timepoint)%>%
#   mutate(Fold = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
#   select(Range, `Subject ID`, Timepoint, Fold) %>% group_by(Timepoint, `Subject ID`) %>%
#   pivot_wider(names_from = Range, values_from = Fold) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "PrototypeOnly_FoldChange_Flow.xlsx"))
# #####

#####
# #plot % igg with all donors
# #look at fold change
# allFlow <- flow %>%
#   mutate(TimeInf = paste0(Timepoint, "_", infect_flag)) %>%
#   group_by(`Subject ID`) %>%
#   mutate(Range = case_when("1_Y" %in% TimeInf ~ "Uh oh",
#                            "15_Y" %in% TimeInf ~ "Uh oh",
#                            "90_Y" %in% TimeInf ~ "Infected Days 15-90",
#                            "180_Y" %in% TimeInf ~ "Infected Days 90-180",
#                            TRUE ~ "Uninfected")) %>%
#   filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
#   mutate(Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180")))
# 
# #cross-reactives
# calcs <- allFlow %>%
#   group_by(Range, Timepoint) %>%
#   summarize(n = n(),
#             mean = mean(ProtoOmi),
#             sd = sd(ProtoOmi)) %>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("% of IgG")+
#   ggtitle("Prototype+BA.1+")+
#   ylim(0, 4)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0.5),
#         legend.position = "none")
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow_percentigg.png"), width = 1.3, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow_percentigg.svg"), width = 1.3, height =1.8)
# 
# allFlow %>% 
#   group_by(Range, `Subject ID`, Timepoint) %>%
#   select(Range, `Subject ID`, Timepoint, ProtoOmi) %>%
#   pivot_wider(names_from = Timepoint, values_from = ProtoOmi) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "CrossReactive_FoldChange_Flow_percentigg.xlsx"))
# 
# ###################BA.1-specific
# calcs <- allFlow %>%
#   group_by(Range, Timepoint) %>%
#   summarize(n = n(),
#             mean = mean(OmiNotProto),
#             sd = sd(OmiNotProto)) %>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("% of IgG")+
#   ggtitle("Prototype-BA.1+")+
#   ylim(0, 4)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0.5),
#         legend.key.size = unit(0.4, "lines"))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow_percentigg.png"), width = 2.5, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow_percentigg.svg"), width = 2.5, height =1.8)
# 
# allFlow %>% 
#   group_by(Range, `Subject ID`, Timepoint) %>%
#   select(Range, `Subject ID`, Timepoint, OmiNotProto) %>%
#   pivot_wider(names_from = Timepoint, values_from = OmiNotProto) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "OmiNotProto_FoldChange_Flow_percentigg.xlsx"))
# 
# ###################proto-specific?
# calcs <- allFlow %>%
#   group_by(Range, Timepoint) %>%
#   summarize(n = n(),
#             mean = mean(ProtoNotOmicron),
#             sd = sd(ProtoNotOmicron)) %>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = mean))+
#   geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
#   geom_line(aes(group = Range, color = Range))+
#   geom_point(shape = 21, aes(fill = Range))+
#   scale_fill_manual(values = rangeColors)+
#   scale_color_manual(values = rangeColors)+
#   ylab("% of IgG")+
#   ggtitle("Prototype+BA.1-")+
#   ylim(0, 4)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         legend.title = element_blank(),
#         plot.title = element_text(hjust = 0.5, vjust = 0.5),
#         legend.key.size = unit(0.4, "lines"))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow_percentigg.png"), width = 2.5, height =1.8, dpi=1000)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow_percentigg.svg"), width = 2.5, height =1.8)
# 
# allFlow %>% 
#   group_by(Range, `Subject ID`) %>%
#   select(Range, `Subject ID`, Timepoint, ProtoNotOmicron) %>% group_by(Timepoint, `Subject ID`) %>%
#   pivot_wider(names_from = Range, values_from = ProtoNotOmicron) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "PrototypeOnly_FoldChange_Flow_percentigg.xlsx"))
# #####

#####
#plot % igg with only citeseq donors
#look at fold change
allFlow <- flow %>%
  mutate(TimeInf = paste0(Timepoint, "_", infect_flag)) %>%
  group_by(`Subject ID`) %>%
  mutate(Range = case_when("1_Y" %in% TimeInf ~ "Uh oh",
                           "15_Y" %in% TimeInf ~ "Uh oh",
                           "90_Y" %in% TimeInf ~ "Infected Days 15-90",
                           "180_Y" %in% TimeInf ~ "Infected Days 90-180",
                           TRUE ~ "Uninfected")) %>%
  filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
  mutate(Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180"))) %>%
  filter((`Subject ID` %in% seuObj$Subject) | (Range %in% c("Infected Days 15-90", "Infected Days 90-180")))


#cross-reactives
calcs <- allFlow %>%
  group_by(Range, Timepoint) %>%
  summarize(n = n(),
            mean = mean(ProtoOmi),
            sd = sd(ProtoOmi)) %>%
  mutate(se = sd / sqrt(n))

ggplot(calcs, aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
  geom_line(aes(group = Range, color = Range))+
  geom_point(shape = 21, aes(fill = Range))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  ylab("% of IgG")+
  ggtitle("Prototype+BA.1+")+
  ylim(0, 4)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow_percentigg_citeseqonly.png"), width = 1.3, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoOmi_InfectionRange_Flow_percentigg_citeseqonly.svg"), width = 1.3, height =1.8)

allFlow %>% 
  group_by(Range, `Subject ID`, Timepoint) %>%
  select(Range, `Subject ID`, Timepoint, ProtoOmi) %>%
  pivot_wider(names_from = Timepoint, values_from = ProtoOmi) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "CrossReactive_FoldChange_Flow_percentigg_citeseqonly.xlsx"))

###################BA.1-specific
calcs <- allFlow %>%
  group_by(Range, Timepoint) %>%
  summarize(n = n(),
            mean = mean(OmiNotProto),
            sd = sd(OmiNotProto)) %>%
  mutate(se = sd / sqrt(n))

ggplot(calcs, aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
  geom_line(aes(group = Range, color = Range))+
  geom_point(shape = 21, aes(fill = Range))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  ylab("% of IgG")+
  ggtitle("Prototype-BA.1+")+
  ylim(0, 4)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        legend.key.size = unit(0.4, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow_percentigg_citeseqonly.png"), width = 2.5, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "OmiNotProto_InfectionRange_Flow_percentigg_citeseqonly.svg"), width = 2.5, height =1.8)

allFlow %>% 
  group_by(Range, `Subject ID`, Timepoint) %>%
  select(Range, `Subject ID`, Timepoint, OmiNotProto) %>%
  pivot_wider(names_from = Timepoint, values_from = OmiNotProto) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "OmiNotProto_FoldChange_Flow_percentigg_citeseqonly.xlsx"))

###################proto-specific?
calcs <- allFlow %>%
  group_by(Range, Timepoint) %>%
  summarize(n = n(),
            mean = mean(ProtoNotOmicron),
            sd = sd(ProtoNotOmicron)) %>%
  mutate(se = sd / sqrt(n))

ggplot(calcs, aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = Range), width = 0.4)+
  geom_line(aes(group = Range, color = Range))+
  geom_point(shape = 21, aes(fill = Range))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  ylab("% of IgG")+
  ggtitle("Prototype+BA.1-")+
  ylim(0, 4)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5),
        legend.key.size = unit(0.4, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow_percentigg_citeseqonly.png"), width = 2.5, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ProtoNotOmi_InfectionRange_Flow_percentigg_citeseqonly.svg"), width = 2.5, height =1.8)

allFlow %>% 
  group_by(Range, `Subject ID`) %>%
  select(Range, `Subject ID`, Timepoint, ProtoNotOmicron) %>% group_by(Timepoint, `Subject ID`) %>%
  pivot_wider(names_from = Range, values_from = ProtoNotOmicron) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "PrototypeOnly_FoldChange_Flow_percentigg_citeseqonly.xlsx"))
#####

#####
#Citeseq data time! What trends do we see in the CITESeq Data?
stats <- df %>% filter(ClusterLabel != "Plasmablast-like", ClusterLabel != "Naive") %>%
  group_by(InfectionRange, Timepoint, ClusterLabel) %>%
  mutate(ClusterLabel = factor(ClusterLabel, levels = c("Acute Activated", "Intermediate", "Resting IgG", "Resting IgA","Atypical", "Plasmablast-like", "Naive"))) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum(linewidth = 0.4)+
  scale_fill_manual(values = shortColors)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  facet_grid(cols = vars(InfectionRange), labeller = label_wrap_gen(width = 15))+
  theme(legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "lines"),
        axis.title.y = element_text(size=7.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7.5,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(size=7.5, face = "bold"),
        panel.spacing = unit(0.4, "lines"))+
  guides(shape = guide_legend(override.aes = list(size=0.3)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "Figure5_AlluvialPlotOfPopulations_infected.png"),width = 3.5, height = 1.9, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "Figure5_AlluvialPlotOfPopulations_infected.svg"),width = 3.5, height = 1.9, units = "in", device = "svg")
dev.off()

#######show changes in cross-reactive activation and BA.1-specific activation over time
#do similar plots to before- what percentage of activated cross-reactive versus prototype-specific cells do we see?
dfAct <- df %>% filter(ClusterLabel != "Naive", ClusterLabel != "Plasmablast-like") %>%
  mutate(ClusterLabel = case_when(ClusterLabel %in% c("Intermediate", "Acute Activated") ~ "Activated"))

calcs <- dfAct %>%
  mutate(ClusterSpec = paste0(ClusterLabel,"_", adj.ProtoOmi)) %>%
  group_by(InfectionRange, Subject, Timepoint, ClusterSpec) %>%
  summarize(n = n()) %>%
  group_by(InfectionRange) %>%
  complete(Subject, Timepoint, ClusterSpec, fill = list(n = 0)) %>%
  group_by(InfectionRange, Subject, Timepoint)%>%
  mutate(Prop = n / sum(n)) %>%
  filter(!is.na(ClusterSpec)) %>%
  group_by(InfectionRange, Timepoint, ClusterSpec) %>%
  summarize(n = n(),
            mean = mean(Prop),
            sd = sd(Prop)) %>%
  mutate(ClusterLabel = str_extract(ClusterSpec, "Activated"),
         adj.ProtoOmi = case_when(str_detect(ClusterSpec, "Proto\\+Omi\\+") ~ "Cross-Reactive",
                                  str_detect(ClusterSpec, "Proto\\-Omi\\+") ~ "BA.1-Specific",
                                  TRUE ~ "Prototype-Specific"),
         se = sd / sqrt(n)) %>%
  filter(!is.na(ClusterLabel))

#cross reactive
calcs %>% filter(adj.ProtoOmi == "Cross-Reactive") %>%
  ggplot(aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = InfectionRange), width = 0.4)+
  geom_line(aes(group = InfectionRange, color = InfectionRange))+
  geom_point(shape = 21, aes(fill = InfectionRange))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"),
                   labels = c("1", "15", "90", "180"))+
  ylab("Proportion of All Cells")+
  ggtitle("Activated Prototype+BA.1+")+
  ylim(0, 1)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 6),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "CrossReactiveActivated_Proportion_Infected.png"), width = 1.4, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "CrossReactiveActivated_Proportion_Infected.svg"), width = 1.4, height =1.8, dpi=1000)


############ba.1 only
calcs %>% filter(adj.ProtoOmi == "BA.1-Specific") %>%
  ggplot(aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = InfectionRange), width = 0.4)+
  geom_line(aes(group = InfectionRange, color = InfectionRange))+
  geom_point(shape = 21, aes(fill = InfectionRange))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"),
                   labels = c("1", "15", "90", "180"))+
  ylab("Proportion of All Cells")+
  ggtitle("Activated Prototype-BA.1+")+
  ylim(0, 1)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 6),
        legend.key.size = unit(0.4, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "BA1Only_ProportionBA1Specific_Infected.png"), width = 2.5, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "BA1Only_ProportionBA1Specific_Infected.svg"), width = 2.5, height =1.8, dpi=1000)

####prototype only
calcs %>% filter(adj.ProtoOmi == "Prototype-Specific") %>%
  ggplot(aes(x = Timepoint, y = mean))+
  geom_errorbar(aes(ymax = se + mean, ymin = mean - se, color = InfectionRange), width = 0.4)+
  geom_line(aes(group = InfectionRange, color = InfectionRange))+
  geom_point(shape = 21, aes(fill = InfectionRange))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"),
                   labels = c("1", "15", "90", "180"))+
  ylab("Proportion of All Cells")+
  ggtitle("Activated Prototype+BA.1-")+
  ylim(0, 1)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 5),
        legend.key.size = unit(0.4, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "PrototypeOnly_citeseq_ProportionPrototypeSpecific_Infected.png"), width = 2.5, height =1.8, dpi=1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "PrototypeOnly_citeseq_ProportionPrototypeSpecific_Infected.svg"), width = 2.5, height =1.8)

###write stats
#write stats
stats <- dfAct %>%
  mutate(ClusterSpec = paste0(ClusterLabel,"_", adj.ProtoOmi)) %>%
  group_by(InfectionRange, Subject, Timepoint, ClusterSpec) %>%
  summarize(n = n()) %>%
  group_by(InfectionRange) %>%
  complete(Subject, Timepoint, ClusterSpec, fill = list(n = 0)) %>%
  group_by(InfectionRange, Subject, Timepoint)%>%
  mutate(Prop = n / sum(n)) %>%
  filter(!is.na(ClusterSpec)) %>%
  mutate(ClusterLabel = str_extract(ClusterSpec, "Activated"),
         adj.ProtoOmi = case_when(str_detect(ClusterSpec, "Proto\\+Omi\\+") ~ "Cross-Reactive",
                                  str_detect(ClusterSpec, "Proto\\-Omi\\+") ~ "BA.1-Specific",
                                  TRUE ~ "Prototype-Specific")) %>%
  filter(!is.na(ClusterLabel))

stats %>%
  filter(adj.ProtoOmi == "Cross-Reactive") %>% select(InfectionRange, Subject, Timepoint, Prop) %>%
  group_by(Subject, Timepoint) %>% pivot_wider(names_from = InfectionRange, values_from = Prop) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "CrossReactive_Activated_Proportion.xlsx"))

stats %>%
  filter(adj.ProtoOmi == "BA.1-Specific") %>% select(InfectionRange, Subject, Timepoint, Prop) %>%
  group_by(Subject, Timepoint) %>% pivot_wider(names_from = InfectionRange, values_from = Prop) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "BA1only_Activated_Proportion.xlsx"))

stats %>%
  filter(adj.ProtoOmi == "Prototype-Specific") %>% select(InfectionRange, Subject, Timepoint, Prop) %>%
  group_by(Subject, Timepoint) %>% pivot_wider(names_from = InfectionRange, values_from = Prop) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "PrototypeOnly_Activated_Proportion.xlsx"))
###

######
#Evolution plots
df$clone <- paste0("s", df$Subject, "_", df$clone_id, "_1")

evolvingInfected <- evolving %>% filter(!is.na(sig), Infection == "Y") %>% mutate(Booster = str_replace(Booster, "Omicron", "BA.1")) %>%
  mutate(InfectionRange = df$InfectionRange[match(clone_id, df$clone)]) %>%
  arrange(sig)

summary <- evolvingInfected %>%
  group_by(adj.ProtoOmi, InfectionRange) %>%
  summarize(n = n()) %>%
  group_by(InfectionRange) %>%
  mutate(order = case_when(adj.ProtoOmi == "Proto+Omi+" ~ 0.0005,
                           adj.ProtoOmi == "Proto+Omi-" ~ 0.0006,
                           adj.ProtoOmi == "Proto-Omi+" ~ 0.0007))

ggplot(evolvingInfected)+
  geom_jitter(aes(fill = adj.ProtoOmi, alpha = sig, shape = adj.ProtoOmi, x = InfectionRange, y= slope), width = 0.3)+
  scale_shape_manual(values = c("Proto+Omi+" = 21, "Proto+Omi-" = 22, "Proto-Omi+"= 24))+
  scale_fill_manual(values = c("Proto+Omi+" = "#FFA630", "Proto+Omi-" =  "#4DA1A9", "Proto-Omi+" = "#D7E8BA"))+
  scale_color_manual(values = c("Proto+Omi+" = "#FFA630", "Proto+Omi-" =  "#4DA1A9", "Proto-Omi+" = "#D7E8BA"))+
  #scale_y_continuous(position = "right")+
  ylab("Slope")+
  coord_cartesian(ylim = c(-4e-04, 4e-04), clip = "off")+
  ggtitle("Lineage Evolution")+
  scale_alpha_discrete(guide = "none")+
  geom_text(data = summary, mapping = aes(label = n, x = InfectionRange, y = I(1), color = adj.ProtoOmi), size = 2, position = position_dodge(width = 1))+
  theme_classic()+
  guides(
         fill = guide_legend(override.aes = list(shape = 21)))+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        axis.title.x = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ClonalEvolution.png"), width = 2.3, height =2.3)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ClonalEvolution.svg"), width = 2.3, height =2.3)
dev.off()
######

######
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

stats <- crossDF %>% filter(Infection == "Y") %>%
  group_by(InfectionRange, Subject, Timepoint, CloneStatusRefined) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  group_by(InfectionRange) %>%
  complete(Subject, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
  mutate(OfficialBooster = crossDF$OfficialBooster[match(Subject, crossDF$Subject)]) %>%
  group_by(InfectionRange, Timepoint, CloneStatusRefined) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum)) %>%
  filter(mean != 0)

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = CloneStatusRefined), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(InfectionRange), labeller = label_wrap_gen(12))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  # scale_fill_manual(values = c("Expanded Pre-Vax" = "#4558A7",
  #                              "Expanded" = "#7598c0",
  #                              "Single Timepoint" = "#EBBB59",
  #                              "Singlet" = "#F2F3F4"))+
  scale_fill_manual(values = c("Expanded Pre-Vax, Post-Infection" = "#600000",
                               "Expanded Pre-Vax" = "#2e5894",
                               "Expanded Pre- and Post-Infection" = "#BF0000",
                               "Expanded Pre-Infection" = "#7392BA",
                               "Expanded Post-Infection" = "#FF0000",
                               "Single Timepoint Pre-Infection" = "#DBE9F3",
                               "Singlet" = "white"))+
  ylab("Proportion")+
  guides(fill = guide_legend(nrow = 4))+
  theme_classic()+
  theme(legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin = margin()),
        panel.spacing = unit(0.35, "lines"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 7),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ClonalRelatednessOverTime_barplots.png"),width = 3.1, height = 3, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "ClonalRelatednessOverTime_barplots.svg"),width = 3, height = 2.8, units = "in")
dev.off()
#######

#######
#make donuts !!
crossDF2 <- crossDF %>% #filter(Infection == "Y") %>% mutate(Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))) #%>% 
  filter(Subject %in% c("5249544848","4950544848"))

plotList <- list()
for(i in unique(crossDF2$Subject)){
  yuh <- crossDF2 %>% filter(Subject == i)
  
  nonSinglets <- unique(yuh$clone_subject_id[duplicated(yuh$clone_subject_id) | duplicated(yuh$clone_subject_id, fromLast=T)])
  yuh$CloneStatus <- ifelse(yuh$clone_subject_id %in% nonSinglets, yuh$clone_subject_id, "Singlet")
  
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
  
  yuh$CloneStatusRefined <- calcs$lab[match(yuh$CloneStatus, calcs$CloneStatus)]
  yuh$CloneStatusRefined <- factor(yuh$CloneStatusRefined, levels = c("Expanded Pre-Vax, Post-Infection",
                                                                      "Expanded Pre- and Post-Infection",
                                                                      "Expanded Post-Infection",
                                                                      "Expanded Pre-Vax",
                                                                      "Expanded Pre-Infection",
                                                                      "Single Timepoint Pre-Infection",
                                                                      "Singlet"))
  
  calcs <- yuh %>%
    group_by(InfectionRange, Subject, Timepoint, CloneStatus) %>%
    summarize(n= n()) %>%
    mutate(Proportion = n / sum(n),
           Total = sum(n),
           ClonalOverlap = yuh$CloneStatusRefined[match(CloneStatus, yuh$CloneStatus)],
           #ClonalOverlap = ifelse(is.na(ClonalOverlap), "Singlet", ClonalOverlap),
           ClonalOverlap = factor(ClonalOverlap, levels = c("Expanded Pre-Vax, Post-Infection",
                                                            "Expanded Pre- and Post-Infection",
                                                            "Expanded Post-Infection",
                                                            "Expanded Pre-Vax",
                                                            "Expanded Pre-Infection",
                                                            "Single Timepoint Pre-Infection",
                                                            "Singlet")))%>%
    arrange(Timepoint, ClonalOverlap)%>%
    mutate(ymax = cumsum(Proportion),
           ymin = c(0, head(ymax, n=-1)))
  
  uniqueMetaVals <- unique(yuh$Timepoint)
  label <- c()
  for(j in c(1:length(uniqueMetaVals))){
    label[j] <- ifelse(length(unique(calcs$Total[calcs$Timepoint == uniqueMetaVals[j]])) < 1, 0, unique(calcs$Total[calcs$Timepoint == uniqueMetaVals[j]]))
  }
  
  dat_text <- data.frame(label = label, Timepoint = uniqueMetaVals)
  #dat_text$StudyWeek <- factor(dat_text$StudyWeek, levels = c(-5, 3, 4, 7.5, 8, 11, 15.5, 16, 19, 24))
  #dat_text <- dat_text %>% filter(label != 0)
  calcs <- calcs %>% mutate(Timepoint = factor(Timepoint, levels = c("Day 0", 'Day 15', "Day 90", "Day 180")))
  dat_text <- dat_text %>% mutate(Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))
  
  p <- ggplot(calcs)+
    geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ClonalOverlap))+
    coord_polar(theta="y")+
    xlim(c(2,4))+
    scale_fill_manual(values = c("Expanded Pre-Vax, Post-Infection" = "#600000",
                                 "Expanded Pre-Vax" = "#2e5894",
                                 "Expanded Pre- and Post-Infection" = "#BF0000",
                                 "Expanded Pre-Infection" = "#7392BA",
                                 "Expanded Post-Infection" = "#FF0000",
                                 "Single Timepoint Pre-Infection" = "#DBE9F3",
                                 "Singlet" = "white"))+
    ggtitle(paste(unique(calcs$InfectionRange), i))+
    guides(fill = "none")+
    geom_text(data = dat_text,
              mapping = aes(x=-Inf, y=-Inf, label = label),
              hjust = 0.5,
              vjust = 0.5,
              size = 2)+
    facet_grid(cols=vars(Timepoint))+
    theme_void()+
    theme(plot.title = element_text(size = 3, hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text(size = 6),
          panel.spacing = unit(-0.49999, "lines"))
  plotList[[i]] <- p
}

#plotList <- plotList[c("4955534848", "5053564848", "4848544848")]
g <- arrangeGrob(grobs = plotList, nrow = length(plotList), ncol=1)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedClonalDonuts.png"), width = 3, height = 1.8, dpi = 1500)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedClonalDonuts.svg"), width = 3, height = 1.8)
#####

#####
#plot clonal lineages over time
clones <- crossDF %>%
  group_by(clone_subject_id) %>%
  mutate(relative = case_when(unique(InfectionRange) == "Infected Days 15-90" & Timepoint %in% c("Day 90", "Day 180") ~ "Post-infection",
                              unique(InfectionRange) == "Infected Days 90-180" & Timepoint %in% c("Day 180") ~ "Post-infection",
                              TRUE ~ "Pre-infection"),
         mu_freq = mu_freq * 100)%>%
  group_by(InfectionRange, clone_subject_id, relative) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  group_by(clone_subject_id) %>%
  mutate(Present = length(unique(relative)) == 2) %>%
  filter(Present)

ggplot(clones, aes(x = relative, y = mean, fill = InfectionRange))+
  geom_line(aes(color = InfectionRange, group = clone_subject_id), alpha= 0.2, linewidth = 0.5)+
  geom_violin()+
  #geom_point(shape = 21)+
  geom_boxplot(aes(fill = InfectionRange),width=0.2, outlier.size = 0)+
  ylab("Mean % VH Mutation")+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  scale_x_discrete(limits = c("Pre-infection", "Post-infection"))+
  facet_grid(cols= vars(InfectionRange))+
  theme_classic()+
  theme(text = element_text(size = 7),
        strip.text = element_text(size = 6),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        axis.title.x = element_blank())
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedSHM_by_clone.png"),width = 1.5, height = 1.8, device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedSHM_by_clone.svg"),width = 1.5, height = 1.8)
dev.off()

#write stats
stats <- crossDF %>%
  group_by(clone_subject_id) %>%
  mutate(relative = case_when(unique(InfectionRange) == "Infected Days 15-90" & Timepoint %in% c("Day 90", "Day 180") ~ "Post-infection",
                              unique(InfectionRange) == "Infected Days 90-180" & Timepoint %in% c("Day 180") ~ "Post-infection",
                              TRUE ~ "Pre-infection"),
         mu_freq = mu_freq * 100)%>%
  group_by(InfectionRange, clone_subject_id, relative) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  group_by(clone_subject_id) %>%
  mutate(Present = length(unique(relative)) == 2) %>%
  filter(Present) %>%
  select(!c(n, Present)) %>% group_by(InfectionRange, clone_subject_id) %>%
  pivot_wider(names_from = relative, values_from = mean) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "PreAndPostInfection_SHMClones.xlsx"))
#####

#####
#look at lineages over time to see better time resolution
clones <- crossDF %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(InfectionRange, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq))

ggplot(clones, aes(x = Timepoint, y = mean, fill = InfectionRange))+
  geom_line(aes(color = InfectionRange, group = clone_subject_id), alpha= 0.5, linewidth = 0.5)+
  geom_point(shape = 21)+
  ylab("Mean % VH Mutation")+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols= vars(InfectionRange))+
  theme_classic()+
  theme(text = element_text(size = 7),
        strip.text = element_text(size = 6),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        axis.title.x = element_blank())
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedSHM_by_clone_alltimes.png"),width = 2.3, height = 1.7, device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "InfectedSHM_by_clone_alltimes.svg"),width = 2.3, height = 1.7)
dev.off()

#write stats
clones <- crossDF %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(InfectionRange, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  select(!n) %>% pivot_wider(names_from = Timepoint, values_from = mean)%>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "SHMClonesOverTime.xlsx"))
#####

#####
#Compare d15 vs d180 for infected vs uninfected donors
compare <- df %>%
          group_by(clone_subject_id) %>%
          mutate(Present = case_when(length(intersect(unique(Timepoint), c("Day 15", "Day 180"))) == 2 ~ "Present",
                                     TRUE ~ "Not Present")) %>%
          filter(Present == "Present", Timepoint %in% c("Day 15", "Day 180"), InfectionRange != "Infected Days 90-180") %>%
          group_by(InfectionRange, clone_subject_id, Timepoint) %>%
          summarize(n = n(),
                    mean = mean(mu_freq * 100))
          
compare %>% mutate(InfectionRange = factor(InfectionRange, levels = c("Uninfected", "Infected Days 15-90"))) %>%
ggplot(aes(x = Timepoint, y = mean, fill = InfectionRange))+
  geom_line(aes(group = clone_subject_id, color = InfectionRange), alpha = 0.8, linewidth = 0.4)+
  geom_violin()+
  geom_boxplot(width = 0.2)+
  ggtitle("Per Lineage")+
  ylab("Mean VH % Mutation")+
  facet_grid(cols = vars(InfectionRange))+
  scale_fill_manual(values = rangeColors)+
  scale_color_manual(values = rangeColors)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.position = "None",
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 5))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "D15_vs_D180_SHM_InfectedvsUninfected.png"), width = 2, height = 1.6, dpi = 1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "D15_vs_D180_SHM_InfectedvsUninfected.svg"), width = 2, height = 1.6, dpi = 1000)

  
#write stats file
df %>%
  group_by(clone_subject_id) %>%
  mutate(Present = case_when(length(intersect(unique(Timepoint), c("Day 15", "Day 180"))) == 2 ~ "Present",
                             TRUE ~ "Not Present")) %>%
  filter(Present == "Present", Timepoint %in% c("Day 15", "Day 180"), InfectionRange != "Infected Days 90-180") %>%
  group_by(InfectionRange, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq * 100)) %>% select(!n) %>%
  group_by(InfectionRange, clone_subject_id) %>%
  pivot_wider(names_from = Timepoint, values_from = mean) %>%
write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "d1vsd15_MeanClonalSHM.xlsx"))
#####


#####
#Plot relative proportions of isotypes over time
#generate colorbrewer palette
# 
# ccall <- c("IGHA1" = "#B2182B", "IGHA2" = "#D6604D",
#            "IGHG1" = "")
# 
# calcs <- df %>% filter(!is.na(c_call)) %>% group_by(InfectionRange, Timepoint, c_call) %>%
#           summarize(n = n()) %>% mutate(Proportion = n/ sum(n))
# 
# ggplot(calcs, aes(x = Timepoint, y = Proportion))+
#   geom_bar(aes(fill = c_call), color = "black", linewidth = 0.1,position = "stack", stat = "identity")+
#   ylab("Proportion")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   scale_
#   theme_classic()+
#   facet_grid(cols = vars(InfectionRange))+
#   theme(text = element_text(size = 6),
#         strip.background = element_blank())
#####
  
#####
#genetics of B cells pre- and post-infection
# crossDF <- df %>% filter(Infection == "Y") %>% filter(adj.ProtoOmi == "Proto+Omi+")
# 
# nonSinglets <- unique(crossDF$clone_subject_id[duplicated(crossDF$clone_subject_id) | duplicated(crossDF$clone_subject_id, fromLast=T)])
# crossDF$CloneStatus <- ifelse(crossDF$clone_subject_id %in% nonSinglets, crossDF$clone_subject_id, "Singlet")
# 
# calcs <- crossDF %>%
#   group_by(CloneStatus, InfectionRange, Timepoint) %>%
#   summarize(n = n()) %>%
#   group_by(CloneStatus, InfectionRange) %>%
#   mutate(relative = case_when(unique(InfectionRange) == "Infected Days 15-90" & Timepoint %in% c("Day 90", "Day 180") ~ "Post-infection",
#                               unique(InfectionRange) == "Infected Days 90-180" & Timepoint %in% c("Day 180") ~ "Post-infection",
#                               TRUE ~ "Pre-infection")) %>%
#   mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
#                          length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) & "Post-infection" %in% unique(relative) ~ "Expanded Pre-Vax, Post-Infection",
#                          length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
#                          length(unique(Timepoint)) > 1 & length(unique(relative)) == 2 ~ "Expanded Pre- and Post-Infection",
#                          length(unique(Timepoint)) > 1 & length(unique(relative)) == 1 & "Pre-infection" %in% unique(relative) ~ "Expanded Pre-Infection",
#                          length(unique(relative)) == 1 & "Post-infection" %in% unique(relative) ~ "Expanded Post-Infection",
#                          length(unique(relative)) == 1 & !("Post-infection" %in% unique(relative)) ~ "Single Timepoint Pre-Infection",
#                          TRUE ~ "CHECK"))
# 
# crossDF$CloneStatusRefined <- calcs$lab[match(crossDF$CloneStatus, calcs$CloneStatus)]
# crossDF$CloneStatusRefined <- factor(crossDF$CloneStatusRefined, levels = c("Expanded Pre-Vax, Post-Infection",
#                                                                             "Expanded Pre- and Post-Infection",
#                                                                             "Expanded Post-Infection",
#                                                                             "Expanded Pre-Vax",
#                                                                             "Expanded Pre-Infection",
#                                                                             "Single Timepoint Pre-Infection",
#                                                                             "Singlet"))
# 
# compareDF <- crossDF %>% filter(CloneStatusRefined %in% c("Expanded Post-Infection", "Expanded Pre-Vax, Post-Infection",
#                                                           "Expanded Pre-Infection", "Expanded Pre- and Post-Infection")) %>%
#   mutate(CloneStatus2 = case_when(CloneStatusRefined == "Expanded Post-Infection" ~ "Post-Infection",
#                                   TRUE ~ "Pre-Infection"))
# 
# #compare SHM overall
# ggplot(compareDF, aes(x = CloneStatus2, y = mu_freq*100, fill = CloneStatus2))+
#   geom_violin()+
#   geom_boxplot(outlier.size = 0, width = 0.1)+
#   scale_x_discrete(limits = c("Pre-Infection", "Post-Infection"))+
#   scale_fill_manual(values = c("Post-Infection" = "#FFAE42",
#                                "Pre-Infection" = "#A8415B"))+
#   ylab("% VH Mutation")+
#   theme_classic()+
#   theme(text = element_text(size = 9, color = "black"),
#         axis.title.x = element_blank(),
#         legend.position = "none",
#         axis.text.x = element_text(angle = 45, hjust =1, vjust = 1, color = "black"),
#         axis.line.x = element_line(linewidth = 0.2),
#         axis.line.y = element_line(linewidth = 0.2),
#         axis.ticks = element_line(linewidth = 0.2))
# ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure 5", "Pre_Post_Infection_SHM.png"), width = 2, height =2.3)
# ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure 5", "Pre_Post_Infection_SHM.svg"), width = 2, height =2.3)
# wilcox.test(compareDF$mu_freq[compareDF$CloneStatus2 == "Post-Infection"],
#             compareDF$mu_freq[compareDF$CloneStatus2 == "Pre-Infection"], paired = FALSE)
# 
# ##########compare vh usage
# vh <- compareDF %>%
#   group_by(Subject, CloneStatus2, v_call) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n)) %>%
#   ungroup()%>%
#   complete(Subject, CloneStatus2, v_call, fill = list(n = 0, Proportion = 0)) %>%
#   group_by(CloneStatus2, v_call) %>%
#   summarize(median = mean(Proportion))
# 
# #plot now
# ggplot(vh, aes(x = v_call, y = CloneStatus2, alpha = median))+
#   geom_tile(color = "white", linewidth = 1, aes(fill = CloneStatus2))+
#   scale_fill_manual(values = c("Post-Infection" = "#FFAE42",
#                                "Pre-Infection" = "#A8415B"), guide = "none")+
#   scale_alpha(limits= c(0,0.16),range = c(0.01, 1))+
#   xlab("VH Gene")+
#   scale_y_discrete(limits = c("Pre-Infection", "Post-Infection"))+
#   theme_classic()+
#   theme(text = element_text(size = 6),
#         legend.key.size = unit(0.3, "lines"),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.y = element_blank(),
#         legend.title = element_blank())
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "VHHeatmap_PreVsPostInfection.png"), width = 7.4, height = 0.9, dpi = 1200)
# #ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S3", "VHHeatmap_CrossvsProto.svg"), width = 7.4, height = 0.9)
# 
# #test statistically?
# vh <- compareDF %>%
#   group_by(Subject, CloneStatus2, v_call) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n)) %>%
#   select(!n) %>%
#   ungroup() %>%
#   complete(Subject, CloneStatus2, v_call, fill = list(Proportion = 0)) %>%
#   pivot_wider(names_from = "CloneStatus2", values_from = "Proportion")
# 
# #do stats and compare VH genes
# uniqueVH <- unique(vh$v_call)
# 
# pvals <- c()
# vhGene <- c()
# for(i in uniqueVH){
#   vhGene <- append(vhGene, i)
#   vhFiltered <- vh %>% filter(v_call == i)
#   
#   pvals <- append(pvals, wilcox.test(vhFiltered$`Pre-Infection`, vhFiltered$`Post-Infection`, paired= TRUE)$p.value)
# }
# 
# pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))
# #####
# 
# #####
# #test vh usage outright between infected/uninfected donors at d180
# d180 <- df %>% filter(adj.ProtoOmi != "Proto+Omi-", Timepoint == "Day 180")
# 
# stats <- d180 %>%
#   group_by(Infection, Subject, v_call) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n)) %>%
#   ungroup()%>%
#   group_by(Infection) %>%
#   complete(Subject, v_call, fill = list(n = 0, Proportion = 0)) %>%
#   group_by(Infection, v_call) %>%
#   summarize(median = mean(Proportion)) %>% mutate(sum = sum(median))
# 
# #plot vh usage
# ggplot(stats, aes(x = v_call, y = Infection, alpha = median))+
#   geom_tile(color = "white", linewidth = 1, aes(fill = Infection))+
#   scale_fill_manual(values = c("Y" = "steelblue",
#                                "N" = "orange"), guide = "none")+
#   scale_alpha(limits= c(0,0.16),range = c(0.01, 1))+
#   xlab("VH Gene")+
#   ylab("Infected?")+
#   scale_y_discrete(limits = c("N", "Y"))+
#   theme_classic()+
#   theme(text = element_text(size = 6),
#         legend.key.size = unit(0.3, "lines"),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         axis.title.y = element_blank(),
#         legend.title = element_blank())
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "VHHeatmap_d180_infected_vs_uninfected.png"), width = 7.4, height = 0.9, dpi = 1200)
# 
# #vh usage testing
# #test statistically?
# vh <- d180 %>%
#   group_by(Subject, v_call) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n)) %>%
#   select(!n) %>%
#   ungroup() %>%
#   complete(Subject, v_call, fill = list(Proportion = 0)) %>%
#   mutate(Infection = df$Infection[match(Subject, df$Subject)])
# 
# #do stats and compare VH genes
# uniqueVH <- unique(vh$v_call)
# 
# pvals <- c()
# vhGene <- c()
# for(i in uniqueVH){
#   vhGene <- append(vhGene, i)
#   vhFiltered <- vh %>% filter(v_call == i)
#   
#   pvals <- append(pvals, wilcox.test(vhFiltered$Proportion[vhFiltered$Infection == "Y"], vhFiltered$Proportion[vhFiltered$Infection == "N"], paired= FALSE)$p.value)
# }
# 
# pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))
# 
# #plot
# vh %>% filter(v_call == "IGHV1-8") %>%
# ggplot(aes(x = Infection, y = Proportion))+
#   geom_point(position = position_jitter(width = 0.2), aes(fill = Infection), shape =21)+
#   scale_fill_manual(values = c("Y" = "steelblue",
#                                "N" = "orange"), guide = "none")+
#   theme_classic()+
#   theme()
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "IGHV1_8_Infection.png"), width = 1.8, height = 2, dpi = 1200)

#####

#####
#Proportion of specificities over time
stats <- df %>%
  group_by(InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         adj.ProtoOmi = factor(adj.ProtoOmi, levels = c("Proto+Omi+",
                                                        "Proto+Omi-","Proto-Omi+"))) %>%
  group_by(InfectionRange) %>%
  complete(Subject, Timepoint, adj.ProtoOmi, fill = list(Proportion = 0))%>%
  group_by(InfectionRange, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n(),
            mean = mean(Proportion),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum),
         InfectionRange = factor(InfectionRange, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180")))
  
ggplot(stats, aes(x = Timepoint, y = mean))+
  geom_bar(stat = "identity", position = "stack", aes(fill = adj.ProtoOmi), color = "black")+
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.3, alpha=0.9) +
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_y_continuous(expand = c(0,0), limits = c(0,1))+
  facet_grid(cols = vars(InfectionRange))+
  ylab("Proportion")+
  theme_classic()+
  theme(text = element_text(size = 9),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 8, angle = 45,hjust = 1, vjust = 1),
        strip.background =element_blank(),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "SpecificityProportionsOverTime.svg"), width = 5, height = 3)

#flow data?
#rerun annotation creation
# allFlow <- flow %>%
#   mutate(TimeInf = paste0(Timepoint, "_", infect_flag)) %>%
#   group_by(`Subject ID`) %>%
#   mutate(Range = case_when("1_Y" %in% TimeInf ~ "Uh oh",
#                            "15_Y" %in% TimeInf ~ "Uh oh",
#                            "90_Y" %in% TimeInf ~ "Infected Days 15-90",
#                            "180_Y" %in% TimeInf ~ "Infected Days 90-180",
#                            TRUE ~ "Uninfected")) %>%
#   filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype/BA.1 mRNA")) %>%
#   mutate(Range = factor(Range, levels = c("Uninfected", "Infected Days 15-90", "Infected Days 90-180"))) %>%
#   filter((`Subject ID` %in% seuObj$Subject) | (Range %in% c("Infected Days 15-90", "Infected Days 90-180")))
# 
# stats <- allFlow %>%
#         select(`Subject ID`, Immunogen, Timepoint, Range, ProtoOmi, ProtoNotOmicron, OmiNotProto) %>%
#         pivot_longer(contains("Proto"), names_to = "Pop", values_to = "PercentIgG") %>%
#         group_by(`Subject ID`, Timepoint, Range) %>%
#         mutate(total = sum(PercentIgG),
#                Prop = PercentIgG / total) %>%
#         group_by(Range, Pop,Timepoint) %>% summarize(Proportion = mean(Prop))
# 
# ggplot(stats, aes(x = Timepoint, y = Proportion))+
#   geom_bar(stat = "identity", position = "stack", aes(fill = Pop), color = "black")+
#   scale_fill_manual(values = c("ProtoOmi" = "#386e72",
#                                "ProtoNotOmicron" = "#95C5C8",
#                                "OmiNotProto" = "#F0C0AA"),
#                     labels = c("ProtoOmi" = "Prototype+BA.1+",
#                                "ProtoNotOmicron" = "Prototype+BA.1-",
#                                "OmiNotProto" = "Prototype-BA.1+"))+
#   facet_grid(cols = vars(Range))+
#   ylab("Proportion")+
#   theme_classic()+
#   theme(text = element_text(size = 9),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 8, angle = 45,hjust = 1, vjust = 1),
#         strip.background =element_blank(),
#         legend.title = element_blank())
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 5", "SpecificityProportionsOverTime_flow.svg"), width = 5, height = 3)

#####