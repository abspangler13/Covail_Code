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

infectionShape = c("Infected" = 21, "Uninfected" = 1)
boostColors = c("BA.1" = "#0F75BD", "Prototype" = "#FBB042")

#read in the data
df <- read_xlsx(here::here("01_raw-data", "MSD_Neut_IgGDataCombined_NoCiteseq.xlsx")) %>%
  filter(!is.na(CorrectedBoost))

###redefine binding
df <- df %>%
  mutate(`Omicron BA.1` = ifelse(`Omicron BA.1.old` < `Omicron BA.1.new`, `Omicron BA.1.old`, `Omicron BA.1.new`)) %>%
  filter(`Omicron BA.1` >= 20000) %>%
  mutate(PrototypeBinding = Prototype >= (`Omicron BA.1.old` * 0.2),
         BA1Binding = `Omicron BA.1` >= 20000,
         BA2Binding = `Omicron BA.2` >= (`Omicron BA.1.new` * 0.2),
         BA45Binding = `Omicron BA 45` >= (`Omicron BA.1.new` * 0.2),
         XBBBinding = `Omicron XBB` >= (`Omicron BA.1.old` * 0.2),
         JN1Binding = `JN1` >= (`Omicron BA.1.old` * 0.2),
         KP2Binding = `KP.2` >= (`Omicron BA.1.new` * 0.2))
#####

#####
#let's compare neut of prototype vs BA.1 --> define "BA.1 bias"
neuts <- df %>%
  mutate(NeutsP = ifelse(str_detect(NeutralizingPopulations, "P\\+"), "Neutralizes", "Doesn't Neutralize"),
         NeutsB = ifelse(str_detect(NeutralizingPopulations, "B\\+"), "Neutralizes", "Doesn't Neutralize")) %>%
  group_by(NeutsP, NeutsB) %>%
  summarize(n = n()) %>%
  mutate(col = case_when(NeutsP == "Neutralizes" & NeutsB == "Neutralizes" ~ "Cross",
                         NeutsP == "Neutralizes" ~ "Prototype",
                         NeutsB == "Neutralizes" ~ "BA.1",
                         TRUE ~ "Neither"))

neuts %>%
  ggplot(aes(y = NeutsP, x = NeutsB))+
  geom_tile(aes(fill = col))+
  xlab("BA.1")+
  ylab("Prototype")+
  geom_text(aes(label = n, color = col), size = 3)+
  scale_fill_manual(values = c("Cross" = "#003049",
                               "Prototype" = "#669BBC",
                               "BA.1" = "#C1121F",
                               "Neither" = "#FDF0D5"))+
  scale_color_manual(values = c("Cross" = "white",
                                "Prototype" = "white",
                                "BA.1" = "white",
                                "Neither" = "black"))+
  ggtitle("Neutralization")+
  theme_classic()+
  theme(text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "BA1BiasedAntibodies.png"),width = 2, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "BA1BiasedAntibodies.svg"),width = 2, height = 2)
######

######
#Compare BA.1-bias to non-bias
neuts <- df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  )) %>%
  mutate(Population = factor(Population, levels = c("Cross-Neutralizing",
                                                    "Prototype-Biased",
                                                    "BA.1-Biased",
                                                    "Non-Neutralizing")))%>%
  group_by(CorrectedBoost, `Infection Status`, donor_id, Population) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>% select(!n)%>%
  group_by(CorrectedBoost, `Infection Status`) %>%
  complete(donor_id, Population, fill = list(Proportion = 0))%>%
  group_by(CorrectedBoost, `Infection Status`, Population) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum),
         BoostInf = paste(CorrectedBoost, ", ", `Infection Status`, sep ="")) %>%
  filter(mean != 0)

ggplot(neuts) +
  geom_bar(aes(x=BoostInf, y=mean, fill = Population), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=BoostInf, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  #facet_grid(cols = vars(InfectionRange), labeller = label_wrap_gen(12))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Prototype, Uninfected", "Prototype, Infected",
                              "BA.1, Uninfected", "BA.1, Infected"))+
  scale_fill_manual(values = c("Cross-Neutralizing" = "#003049",
                               "Prototype-Biased" = "#669BBC",
                               "BA.1-Biased" = "#C1121F",
                               "Non-Neutralizing" = "#FDF0D5"))+
  ylab("Proportion")+
  guides(fill = guide_legend(nrow = 4))+
  theme_classic()+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.key.size = unit(0.3, 'cm'),
        text = element_text(size = 5),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin = margin()),
        panel.spacing = unit(0.35, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "Barplot_NeutralizingPopulations.png"), width = 2.3, height = 2.5)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "Barplot_NeutralizingPopulations.svg"), width = 2.3, height = 2.5)

#plot out proportion of each donor that fits in each population
neuts <- df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  )) %>%
  mutate(Population = factor(Population, levels = c("Cross-Neutralizing",
                                                    "Prototype-Biased",
                                                    "BA.1-Biased",
                                                    "Non-Neutralizing")),
         CorrectedBoost = factor(CorrectedBoost, levels = c("Prototype", "BA.1")),
         `Infection Status` = factor(`Infection Status`, levels = c("Uninfected", "Infected")))%>%
  group_by(CorrectedBoost, `Infection Status`, donor_id, Population) %>%
  summarize(n = n()) %>%
  group_by(CorrectedBoost, `Infection Status`, donor_id) %>%
  complete(Population, fill = list(n = 0))%>%
  mutate(Proportion = n / sum(n),
         BoostInf = paste0(CorrectedBoost, ", ", `Infection Status`))

neuts %>%
  ggplot(aes(x = Population, y = Proportion, color = CorrectedBoost, fill = CorrectedBoost))+
  geom_point(aes(shape = `Infection Status`) ,position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.01), size = 1)+
  ylab("Proportion of mAbs")+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  scale_shape_manual(values = infectionShape)+
  ylim(0,1)+
  theme_classic()+
  guides(shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 5),
        legend.key.size = unit(0.4, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "ProportionOfPopulations_perdonor.png"), width = 2.6, height = 2, dpi = 1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "ProportionOfPopulations_perdonor.svg"), width = 2.6, height = 2)


#write stats sheet
df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  )) %>%
  mutate(Population = factor(Population, levels = c("Cross-Neutralizing",
                                                    "Prototype-Biased",
                                                    "BA.1-Biased",
                                                    "Non-Neutralizing")))%>%
  group_by(CorrectedBoost, `Infection Status`, donor_id, Population) %>%
  summarize(n = n()) %>%
  group_by(CorrectedBoost, `Infection Status`, donor_id) %>%
  complete(Population, fill = list(n = 0))%>%
  mutate(Proportion = n / sum(n),
         BoostInf = paste0(CorrectedBoost, ", ", `Infection Status`)) %>%
  ungroup()%>%
  select(BoostInf, donor_id, Population, Proportion) %>%
  group_by(donor_id, Population) %>%
  pivot_wider(names_from = BoostInf, values_from = Proportion)%>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures","Figure 6", "all_cells", "ProportionOfEachNeutPopulation_allcells.xlsx"))

######

######
#of each category of antibody, do we see improved neut of XBB or JN.1?
neuts <- df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  ),
  JX = case_when(
    str_detect(NeutralizingPopulations, "X\\+J\\+") ~ "Neutralizes Both",
    str_detect(NeutralizingPopulations, "X\\+J-") ~ "Neutralizes XBB",
    str_detect(NeutralizingPopulations, "X-J\\+") ~ "Neutralizes JN.1",
    TRUE ~ "Neutralizes Neither"
  )) %>%
  group_by(Population, JX) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         JX = factor(JX, levels = c(
           "Neutralizes Both",
           "Neutralizes JN.1",
           "Neutralizes XBB",
           "Neutralizes Neither"
         )))

ggplot(neuts, aes(x = Population, y = Proportion, fill = JX))+
  geom_bar(stat = "identity", position = "stack", color= "black")+
  scale_fill_manual(values = c("Neutralizes Both" = "black",
                               "Neutralizes JN.1" = "gray33",
                               "Neutralizes XBB" = "gray66",
                               "Neutralizes Neither" = "whitesmoke"))+
  ylab("Proportion")+
  ggtitle("XBB.1.5 and JN.1 Neutralization")+
  scale_x_discrete(limits = c("Cross-Neutralizing", "BA.1-Biased", "Prototype-Biased", "Non-Neutralizing"))+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "NeutralizingBias_BreadthTowardsXBBJN1.png"), width = 2.8, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "NeutralizingBias_BreadthTowardsXBBJN1.svg"), width = 2.8, height = 2)

#write numbers for Sarah
write_xlsx(neuts,here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "XBBJN1Neut_ByBias.xlsx"))

#also do by donor
df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  ),
  JX = case_when(
    str_detect(NeutralizingPopulations, "X\\+J\\+") ~ "Neutralizes Both",
    str_detect(NeutralizingPopulations, "X\\+J-") ~ "Neutralizes XBB",
    str_detect(NeutralizingPopulations, "X-J\\+") ~ "Neutralizes JN.1",
    TRUE ~ "Neutralizes Neither"
  )) %>%
  group_by(CorrectedBoost, `Infection Status`,`Donor ID`,Population, JX) %>%
  summarize(n = n()) %>%
  mutate(JX = factor(JX, levels = c(
           "Neutralizes Both",
           "Neutralizes JN.1",
           "Neutralizes XBB",
           "Neutralizes Neither"
         ))) %>%
  pivot_wider(names_from = JX, values_from=n)%>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "XBBJN1Neut_ByBias_perdonor.xlsx"))
######




#####
#What if we combine the infecteds and see what happens
neuts <- df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  )) %>%
  mutate(Population = factor(Population, levels = c("Cross-Neutralizing",
                                                    "Prototype-Biased",
                                                    "BA.1-Biased",
                                                    "Non-Neutralizing")),
         CorrectedBoost = factor(CorrectedBoost, levels = c("Prototype", "BA.1")),
         `Infection Status` = factor(`Infection Status`, levels = c("Uninfected", "Infected")))%>%
  group_by(CorrectedBoost, donor_id, Population) %>%
  summarize(n = n()) %>%
  group_by(CorrectedBoost, donor_id) %>%
  complete(Population, fill = list(n = 0))%>%
  mutate(Proportion = n / sum(n))

neuts %>%
  ggplot(aes(x = Population, y = Proportion, color = CorrectedBoost, fill = CorrectedBoost))+
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.01), size = 1)+
  ylab("Proportion of mAbs")+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  #scale_shape_manual(values = infectionShape)+
  ylim(0,1)+
  theme_classic()+
  guides(shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 5),
        legend.key.size = unit(0.4, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank())
#ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "ProportionOfPopulations_perdonor.png"), width = 2.6, height = 2, dpi = 1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "all_cells", "ProportionOfPopulations_perdonor_combinedbooster.svg"), width = 2.6, height = 2)

#calcs
neuts <- df %>%
  mutate(Population = case_when(
    str_detect(NeutralizingPopulations, "P\\+B\\+") ~ "Cross-Neutralizing",
    str_detect(NeutralizingPopulations, "P\\+B-") ~ "Prototype-Biased",
    str_detect(NeutralizingPopulations, "P-B\\+") ~ "BA.1-Biased",
    TRUE ~ "Non-Neutralizing"
  )) %>%
  mutate(Population = factor(Population, levels = c("Cross-Neutralizing",
                                                    "Prototype-Biased",
                                                    "BA.1-Biased",
                                                    "Non-Neutralizing")),
         CorrectedBoost = factor(CorrectedBoost, levels = c("Prototype", "BA.1")),
         `Infection Status` = factor(`Infection Status`, levels = c("Uninfected", "Infected")))%>%
  group_by(CorrectedBoost, donor_id, Population) %>%
  summarize(n = n()) %>%
  group_by(CorrectedBoost, donor_id) %>%
  complete(Population, fill = list(n = 0))%>%
  mutate(Proportion = n / sum(n))

#cross neutralizing
wilcox.test(neuts$Proportion[neuts$CorrectedBoost== "Prototype" & neuts$Population == "Cross-Neutralizing"], neuts$Proportion[neuts$CorrectedBoost == "BA.1"& neuts$Population == "Cross-Neutralizing"], paired = FALSE)

#prototype biased
wilcox.test(neuts$Proportion[neuts$CorrectedBoost== "Prototype" & neuts$Population == "Prototype-Biased"], neuts$Proportion[neuts$CorrectedBoost == "BA.1"& neuts$Population == "Prototype-Biased"], paired = FALSE)

#BA.1 biased
wilcox.test(neuts$Proportion[neuts$CorrectedBoost== "Prototype" & neuts$Population == "BA.1-Biased"], neuts$Proportion[neuts$CorrectedBoost == "BA.1"& neuts$Population == "BA.1-Biased"], paired = FALSE)

#Non-neutralizing
wilcox.test(neuts$Proportion[neuts$CorrectedBoost== "Prototype" & neuts$Population == "Non-Neutralizing"], neuts$Proportion[neuts$CorrectedBoost == "BA.1"& neuts$Population == "Non-Neutralizing"], paired = FALSE)

#####