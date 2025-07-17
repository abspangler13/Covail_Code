library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(RColorBrewer)

#load in the flow data and make the necessary variables
flow <- read.csv(here::here("04_Analysis", "data_objects", "figure_testing", "ratpig_phenotyping", "sort1_phenotyping.wsp FlowJo table.csv")) %>%
        filter(!is.na(Activated_Atypical)) %>% select(!X.1) %>%
        mutate(ID = str_extract(X, "(?<=Specimen_001_)[0-9]+"))

infected <- c("206230717", "206231302", "206311727", "206333104", "206362879","206363191","206382618","206383247","206383617","206447647")
flow$Infected <- flow$ID %in% infected

#########
#Rough acute activated cluster
######activated
ggplot(flow, aes(x = Infected, y = Activated_Atypical, shape = Infected))+
  geom_jitter(width = 0.1, fill = "black")+
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 21))+
  scale_x_discrete(labels = c("TRUE" = "Infected", "FALSE" = "Uninfected"))+
  theme_classic()+
  ylab("Percentage of IgG+")+
  ylim(0, 17)+
  theme(text = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "figure_testing", "ratpig_phenotyping", "activated_atypical.svg"), width = 2, height = 3)
wilcox.test(flow$Activated_Atypical[flow$Infected], flow$Activated_Atypical[!flow$Infected], paired = FALSE)

#####antigen-specific activated
ggplot(flow, aes(x = Infected, y = AntigenSpecific_ActivatedAtypical, shape = Infected))+
  geom_jitter(width = 0.1, fill = "black")+
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 21))+
  scale_x_discrete(labels = c("TRUE" = "Infected", "FALSE" = "Uninfected"))+
  theme_classic()+
  ylab("Percentage of BA.1++")+
  #ylim(0, 17)+
  theme(text = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "figure_testing", "ratpig_phenotyping", "activated_atypical_antigenspecific.svg"), width = 2, height = 3)
wilcox.test(flow$AntigenSpecific_ActivatedAtypical[flow$Infected], flow$AntigenSpecific_ActivatedAtypical[!flow$Infected], paired = FALSE)

########freq of igg???
ggplot(flow, aes(x = Infected, y = AntigenSpecificActivated_freqofigg, shape = Infected))+
  geom_jitter(width = 0.1, fill = "black")+
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 21))+
  scale_x_discrete(labels = c("TRUE" = "Infected", "FALSE" = "Uninfected"))+
  theme_classic()+
  ylab("Percentage of IgG+")+
  #ylim(0, 17)+
  theme(text = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "figure_testing", "ratpig_phenotyping", "activated_atypical_antigenspecific_freqofigg.svg"), width = 2, height = 3)
wilcox.test(flow$AntigenSpecificActivated_freqofigg[flow$Infected], flow$AntigenSpecificActivated_freqofigg[!flow$Infected], paired = FALSE)
######

######
#Acute activated cluster- factoring in CD27+
########freq of igg???
ggplot(flow, aes(x = Infected, y = Activated_CD27Pos, shape = Infected))+
  geom_jitter(width = 0.1, fill = "black")+
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 21))+
  scale_x_discrete(labels = c("TRUE" = "Infected", "FALSE" = "Uninfected"))+
  theme_classic()+
  ylab("Percentage of IgG+")+
  ggtitle("CD21loCD27+")+
  ylim(0, 5)+
  theme(text = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "figure_testing", "ratpig_phenotyping", "activated_notatypical.svg"), width = 2, height = 3)
wilcox.test(flow$Activated_CD27Pos[flow$Infected], flow$Activated_CD27Pos[!flow$Infected], paired = FALSE)
#####

######
#intermediate?
########
ggplot(flow, aes(x = Infected, y = CD21hiCD11cpos, shape = Infected))+
  geom_jitter(width = 0.1, fill = "black")+
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 21))+
  scale_x_discrete(labels = c("TRUE" = "Infected", "FALSE" = "Uninfected"))+
  theme_classic()+
  ylab("Percentage of IgG+")+
  ggtitle("CD21loCD27+")+
  ylim(0, 5)+
  theme(text = element_text(color = "black", size = 8),
        axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 7, hjust = 0.5),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "figure_testing", "ratpig_phenotyping", "activated_notatypical.svg"), width = 2, height = 3)
wilcox.test(flow$Activated_CD27Pos[flow$Infected], flow$Activated_CD27Pos[!flow$Infected], paired = FALSE)
#####