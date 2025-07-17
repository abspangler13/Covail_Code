library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(shazam)

#####
#set colors
#let's work on this
immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#068041",
                     "Prototype + Beta" = "#c47002",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f",
                     "Delta + BA.1" = "firebrick4",
                     "Beta + BA.1" = "wheat1")

infectionShape = c("Infected" = 21, "Uninfected" = 1)
boostColors = c("BA.1" = "#2AB673", "Prototype" = "#FBB042")

#read in the data
totaldf <- read_xlsx(here::here("01_raw-data", "MSD_Neut_IgGDataCombined_NoCiteseq.xlsx")) %>%
  filter(!is.na(CorrectedBoost)) #%>%
# mutate(CorrectedBoost = factor(CorrectedBoost, levels = c("Prototype", "BA.1")),
#        `Infection Status` = factor(`Infection Status`, levels = c("Uninfected", "Infected")))

###redefine binding
df <- totaldf %>%
  mutate(`Omicron BA.1` = ifelse(`Omicron BA.1.old` < `Omicron BA.1.new`, `Omicron BA.1.old`, `Omicron BA.1.new`)) %>%
  filter(`Omicron BA.1` >= 20000) %>%
  mutate(PrototypeBinding = Prototype >= (`Omicron BA.1.old` * 0.2),
         BA1Binding = `Omicron BA.1` >= 20000,
         BA2Binding = `Omicron BA.2` >= (`Omicron BA.1.new` * 0.2),
         BA45Binding = `Omicron BA 45` >= (`Omicron BA.1.new` * 0.2),
         XBBBinding = `Omicron XBB` >= (`Omicron BA.1.old` * 0.2),
         JN1Binding = `JN1` >= (`Omicron BA.1.old` * 0.2),
         KP2Binding = `KP.2` >= (`Omicron BA.1.new` * 0.2),
         `Infection Status` = factor(`Infection Status`, levels = c("Uninfected", "Infected")),
         CorrectedBoost = factor(CorrectedBoost, levels = c("Prototype", "BA.1")),
         `Donor ID` = donor_id)

#####import raw dataset
raw <- read_xlsx(here::here("01_raw-data", "COVAIL_unfiltered_MSD_Neut_IgG_Data.xlsx"))

#####

#####
#raw plots 
rawcalcs <- raw %>%
  select(ID, IgGQuant, Prototype, `Omicron BA.1`, `Omicron BA.2`, `Omicron BA 45`, `Omicron XBB`, JN1, KP.2) %>%
  pivot_longer(cols = c(contains("Omicron"), Prototype, JN1, KP.2), names_to = "RBD") %>%
  mutate(RBD = str_remove(RBD, "Omicron "),
         RBD = case_when(RBD == "BA 45" ~ "BA.4/5",
                         RBD == "JN1" ~ "JN.1",
                         RBD == "XBB" ~ "XBB.1.5",
                         TRUE ~ RBD),
         RBD = factor(RBD, levels = c("Prototype","BA.1", "BA.2", "BA.4/5", "XBB.1.5", "JN.1", "KP.2"))) %>% filter(RBD == "BA.1")

ggplot(rawcalcs, aes(x = value))+
  geom_histogram(fill = "white", color = "black", binwidth = 50000, linewidth = 0.2)+
  geom_vline(xintercept = 20000, color = "red", linetype = 2)+
  xlab("ECL Units")+
  #scale_x_continuous(trans=scales::pseudo_log_trans(sigma = 1))+
  #facet_grid(rows = vars(RBD))+
  ggtitle("BA.1 Binding")+
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust =1),
        axis.line = element_line(linewidth = 0.2),
        axis.title.y = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "BA1_ECLUDistribution.svg"), width = 4, height = 3)


######
#Now show distributions of all but BA1 post-processing
calcs <- df %>%
  select(CELL, IgGQuant, Prototype, `Omicron BA.1`, `Omicron BA.2`, `Omicron BA 45`, `Omicron XBB`, JN1, KP.2) %>%
  pivot_longer(cols = c(contains("Omicron"), Prototype, JN1, KP.2), names_to = "RBD") %>%
  mutate(RBD = str_remove(RBD, "Omicron "),
         RBD = case_when(RBD == "BA 45" ~ "BA.4/5",
                         RBD == "JN1" ~ "JN.1",
                         RBD == "XBB" ~ "XBB.1.5",
                         TRUE ~ RBD),
         RBD = factor(RBD, levels = c("Prototype","BA.1", "BA.2", "BA.4/5", "XBB.1.5", "JN.1", "KP.2"))) %>% filter(RBD != "BA.1")

ggplot(calcs, aes(x = value))+
  geom_histogram(fill = "white", color = "black", binwidth = 50000, linewidth = 0.2)+
  xlab("ECL Units")+
  #scale_x_continuous(trans=scales::pseudo_log_trans(sigma = 1))+
  scale_y_continuous(expand = c(0, 0))+
  facet_grid(rows = vars(RBD))+
  geom_hline(yintercept = 0)+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1, vjust =1),
        axis.line = element_line(linewidth = 0.2),
        axis.title.y = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AllRBDs_ECLUDistribution_filtered_data.svg"), width = 4, height = 6)
