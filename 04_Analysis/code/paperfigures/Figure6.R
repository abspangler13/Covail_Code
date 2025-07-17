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

#df with smartseq
citeRat <- read_xlsx(here::here("01_raw-data", "CombinedMSD_Neut_Smartseq_OverlappingCITESeq.xlsx"))
#####

#####
#Make a barplot showing imprinting
stats <- df %>%
  group_by(CorrectedBoost, `Infection Status`, `Donor ID`, PrototypeBinding) %>%
  summarize(n = n()) %>%
  group_by(CorrectedBoost, `Infection Status`) %>%
  complete(`Donor ID`, PrototypeBinding, fill = list(Proportion = 0, n=0)) %>%
  mutate(Status = paste0(CorrectedBoost, ", ", `Infection Status`),
         Status = factor(Status, levels = c("Prototype, Uninfected", "Prototype, Infected",
                                            "BA.1, Uninfected", "BA.1, Infected")),
         donor_id = as.character(`Donor ID`)) %>% rename("Prototype Binding"="PrototypeBinding")

ggplot(stats, aes(x = donor_id, y = n))+
  geom_bar(stat = "identity",position = "stack", aes(fill = `Prototype Binding`), color = "black", linewidth = 0.1)+
  #facet_grid(cols = vars(`Infection Status`), rows = vars(CorrectedBoost), scale = "free", axes = "all_x")+
  facet_wrap(vars(Status), scale = "free_x")+
  ylab("Number of mAbs")+
  xlab("Donor")+
  scale_fill_manual(values = c("FALSE" = "#C5A586", "TRUE" = "#EEE9E3"), labels = c("TRUE" = "Prototype+BA.1+", "FALSE" = "Prototype-BA.1+"))+
  theme_classic()+
  theme(text = element_text(size = 12, color = "black"),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.6, "lines"))
ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure 6", "PerDonor_DeNovoResponse.png"), width = 4.8, height = 4.5)
ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure 6", "PerDonor_DeNovoResponse.svg"), width = 4.8, height = 4.5)

#create a table
df %>%
  group_by(`Infection Status`, CorrectedBoost, `Donor ID`, PrototypeBinding) %>%
  summarize(n = n()) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "Subject_PrototypeBinding_table.xlsx"))
#####

#####
#VH mutation of de novo response vs recall
citeRat <- citeRat %>% filter(!is.na(IgGQuant))

ratpigSHM <- observedMutations(citeRat, sequenceColumn = "sequence_alignment",
                               germlineColumn = "germline_alignment",
                               regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)

ratpigSHM2 <- ratpigSHM %>% mutate(DeNovo = BindingPopulation %in% c("BA.1-specific", "P-B+X+J+", "P-B+X+J-", "P-B+X-J+", "P-B+X-J-")) %>%
  filter(!(CorrectedBoost == "Prototype" & `Infection Status` == "Uninfected" & DeNovo))

ggplot(ratpigSHM2, aes(x= DeNovo, y= mu_freq * 100))+
  geom_violin(aes(fill= DeNovo))+
  geom_boxplot(aes(fill = DeNovo), width = 0.2)+
  ylab("% VH Mutation")+
  scale_fill_manual(values = c("TRUE" = "#C5A586", "FALSE" = "#EEE9E3"))+
  scale_x_discrete(limits = c("FALSE", "TRUE"), labels = c("TRUE" = "Prototype-BA.1+", "FALSE" = "Prototype+BA.1+"))+
  theme_classic()+
  theme(text = element_text(size = 10, color = "black"),
        axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "Denovo_Vs_Recall_SHM.png"), width = 2.2, height = 3)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "Denovo_Vs_Recall_SHM.svg"), width = 2.2, height = 3)

denovo <- ratpigSHM2$mu_freq[ratpigSHM2$DeNovo]
recall <- ratpigSHM2$mu_freq[!ratpigSHM2$DeNovo]

wilcox.test(denovo, recall, paired= FALSE)
#####

#####
#Binding and neut graphs now
###################################By infection
#show neuts vs binding from Kyle's graph
data <- df %>% mutate(
  NeutsPro = stringr::str_count(NeutralizingPopulations, "P\\+"),
  NeutsBA1 = stringr::str_count(NeutralizingPopulations, "B\\+"),
  NeutsJN1 = stringr::str_count(NeutralizingPopulations, "J\\+"),
  NeutsXBB = stringr::str_count(NeutralizingPopulations, "X\\+")
)

#compare binding
bindingdf <- data  |>
  select(`Donor ID`, CorrectedBoost, `Infection Status`, `RATPIg Well`, PrototypeBinding, BA1Binding, BA2Binding, BA45Binding, XBBBinding, JN1Binding, KP2Binding) %>%
  pivot_longer(cols = contains("Binding"), names_to = "Antigen", values_to = "Yes?") %>%
  group_by(CorrectedBoost, `Infection Status`, `Donor ID`, Antigen, `Yes?`) %>%
  summarize(n =n ()) %>%
  mutate(proportion = n / sum(n)) %>%
  filter(`Yes?`) %>%
  mutate(Antigens = factor(Antigen, levels = c("PrototypeBinding", "BA1Binding", "BA2Binding", "BA45Binding", "XBBBinding", "JN1Binding", "KP2Binding"))) 

ggplot(bindingdf, aes(x = factor(Antigen), y = proportion, shape = `Infection Status`))+
  geom_point(aes(fill = CorrectedBoost, group = `Infection Status`, color =CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 1.5, stroke = 0.1)+
  ylim(0,1)+
  theme_classic()+
  #scale_fill_manual(values = c("darkslateblue", "dodgerblue"))+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  ggtitle("Binding")+
  scale_x_discrete("Antigen", labels = c("PrototypeBinding" = "Prototype",
                                         "BA1Binding" = "BA.1",
                                         "BA2Binding" = "BA.2",
                                         "BA45Binding" = "BA.4/5",
                                         "XBBBinding" = "XBB.1.5",
                                         "JN1Binding" = "JN.1",
                                         "KP2Binding" = "KP.2"),
                   limits = c("PrototypeBinding", "BA1Binding", "BA2Binding", "BA45Binding", "XBBBinding",
                              "JN1Binding", "KP2Binding"))+
  scale_shape_manual(values = c("Uninfected" = 21, "Infected" = 1))+
  ylab("Proportion")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "BindingComparison_byInfection.png"), dpi = 1000, width = 2.8, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "BindingComparison_byInfection.svg"), width = 2.8, height = 2)

#write excel file
write_xlsx(bindingdf, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "BindingComparison.xlsx"))

#now neutralization
neutsdf <- data |>
  select(`Donor ID`, CorrectedBoost, `Infection Status`, `RATPIg Well`, NeutsPro, NeutsBA1, NeutsXBB, NeutsJN1) |>
  pivot_longer(cols = contains("Neuts"), names_to = "Antigen", values_to = "Yes?") |>
  group_by(CorrectedBoost, `Infection Status`, `Donor ID`, Antigen, `Yes?`) |>
  summarize(n =n ()) |>
  mutate(proportion = n / sum(n)) |>
  filter(`Yes?` == 1) |>
  mutate(Antigens = case_when(Antigen == "NeutsPro" ~ "Prototype",
                              Antigen == "NeutsBA1" ~ "BA.1",
                              Antigen == "NeutsXBB" ~ "XBB.1.5",
                              Antigen == "NeutsJN1" ~ "JN.1"),
         Antigens = factor(Antigens, levels = c("Prototype", "BA.1", "XBB.1.5", "JN.1"))) |> #factor in mutate so you don't need to specify that dataframe
  group_by(CorrectedBoost, `Infection Status`) |> #
  complete(`Donor ID`, Antigens, fill = list(proportion = 0))

ggplot(neutsdf, aes(x = factor(Antigens), y = proportion, shape = `Infection Status`))+
  geom_point(aes(fill = CorrectedBoost, group = `Infection Status`, color = CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 1.5, stroke = 0.1)+
  ylim(0,1)+
  theme_classic()+
  #scale_fill_manual(values = c("darkslateblue", "dodgerblue"))+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  ggtitle("Neutralization")+
  # scale_x_discrete("Antigens", labels = c("NeutsPro" = "Prototype",
  #                                        "NeutsBA1" = "BA.1",
  #                                        "NeutsXBB" = "XBB.1.5",
  #                                        "NeutsJN1" = "JN.1"),
  #                  limits = c("NeutsPro", "NeutsBA1", "NeutsXBB", "NeutsJN1"))+
  scale_shape_manual(values = c("Uninfected" = 21, "Infected" = 1))+
  ylab("Proportion")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutsComparison_ByInfection.png"), dpi = 1000, width = 2.4, height = 1.88)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutsComparison_ByInfection.svg"), dpi = 1000, width = 2.4, height = 1.88)

#write xlsx file
write_xlsx(neutsdf, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "NeutsComparison.xlsx"))


###################################By vaccination
#compare binding
bindingdf <- data  |>
  select(`Donor ID`, CorrectedBoost, `Infection Status`, `RATPIg Well`, PrototypeBinding, BA1Binding, BA2Binding, BA45Binding, XBBBinding, JN1Binding, KP2Binding) %>%
  pivot_longer(cols = contains("Binding"), names_to = "Antigen", values_to = "Yes?") %>%
  group_by(CorrectedBoost, `Infection Status`, `Donor ID`, Antigen, `Yes?`) %>%
  summarize(n =n ()) %>%
  mutate(proportion = n / sum(n)) %>%
  filter(`Yes?`) %>%
  mutate(Antigens = factor(Antigen, levels = c("PrototypeBinding", "BA1Binding", "BA2Binding", "BA45Binding", "XBBBinding", "JN1Binding", "KP2Binding"))) 

ggplot(bindingdf, aes(x = factor(Antigen), y = proportion, shape = `Infection Status`))+
  geom_point(aes(fill = CorrectedBoost, group = CorrectedBoost, color = CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 1.5, stroke = 0.1)+
  ylim(0,1)+
  theme_classic()+
  #scale_fill_manual(values = c("darkslateblue", "dodgerblue"))+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  ggtitle("Binding")+
  scale_x_discrete("Antigen", labels = c("PrototypeBinding" = "Prototype",
                                         "BA1Binding" = "BA.1",
                                         "BA2Binding" = "BA.2",
                                         "BA45Binding" = "BA.4/5",
                                         "XBBBinding" = "XBB.1.5",
                                         "JN1Binding" = "JN.1",
                                         "KP2Binding" = "KP.2"),
                   limits = c("PrototypeBinding", "BA1Binding", "BA2Binding", "BA45Binding", "XBBBinding",
                              "JN1Binding", "KP2Binding"))+
  scale_shape_manual(values = c("Uninfected" = 21, "Infected" = 1))+
  ylab("Proportion")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 7),
        legend.key.size = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank())
#legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "BindingComparison_byVaccination.png"), dpi = 1000, width = 2.8, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "BindingComparison_byVaccination.svg"), dpi = 1000, width = 2.8, height = 2)

#now neutralization
neutsdf <- data |>
  select(`Donor ID`, CorrectedBoost, `Infection Status`, `RATPIg Well`, NeutsPro, NeutsBA1, NeutsXBB, NeutsJN1) |>
  pivot_longer(cols = contains("Neuts"), names_to = "Antigen", values_to = "Yes?") |>
  group_by(CorrectedBoost, `Infection Status`, `Donor ID`, Antigen, `Yes?`) |>
  summarize(n =n ()) |>
  mutate(proportion = n / sum(n)) |>
  filter(`Yes?` == 1) |>
  mutate(Antigens = case_when(Antigen == "NeutsPro" ~ "Prototype",
                              Antigen == "NeutsBA1" ~ "BA.1",
                              Antigen == "NeutsXBB" ~ "XBB.1.5",
                              Antigen == "NeutsJN1" ~ "JN.1"),
         Antigens = factor(Antigens, levels = c("Prototype", "BA.1", "XBB.1.5", "JN.1"))) |> #factor in mutate so you don't need to specify that dataframe
  group_by(CorrectedBoost, `Infection Status`) |> #
  complete(`Donor ID`, Antigens, fill = list(proportion = 0))

ggplot(neutsdf, aes(x = factor(Antigens), y = proportion, shape = `Infection Status`))+
  geom_point(aes(fill = CorrectedBoost, group = CorrectedBoost, color = CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 1.5, stroke = 0.1)+
  ylim(0,1)+
  theme_classic()+
  #scale_fill_manual(values = c("darkslateblue", "dodgerblue"))+
  scale_fill_manual(values = boostColors)+
  scale_color_manual(values = boostColors)+
  ggtitle("Neutralization")+
  # scale_x_discrete("Antigens", labels = c("NeutsPro" = "Prototype",
  #                                        "NeutsBA1" = "BA.1",
  #                                        "NeutsXBB" = "XBB.1.5",
  #                                        "NeutsJN1" = "JN.1"),
  #                  limits = c("NeutsPro", "NeutsBA1", "NeutsXBB", "NeutsJN1"))+
  scale_shape_manual(values = c("Uninfected" = 21, "Infected" = 1))+
  ylab("Proportion")+
  guides(fill = guide_legend(override.aes = list(shape=21)),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutsComparison_ByVaccination.png"), dpi = 1000, width = 1.6, height = 1.88)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutsComparison_ByVaccination.svg"), dpi = 1000, width = 1.6, height = 1.88)
#######

#####
#Schism the neut-biased populations
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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "Barplot_NeutralizingPopulations.png"), width = 2.3, height = 2.5)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "Barplot_NeutralizingPopulations.svg"), width = 2.3, height = 2.5)

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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "ProportionOfPopulations_perdonor.png"), width = 2.6, height = 2, dpi = 1000)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "ProportionOfPopulations_perdonor.svg"), width = 2.6, height = 2)


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
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures","Figure 6", "ProportionOfEachNeutPopulation_allcells.xlsx"))
########

######
#What proportion of each population neuts XBB and JN1?
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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutralizingBias_BreadthTowardsXBBJN1.png"), width = 2.8, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 6", "NeutralizingBias_BreadthTowardsXBBJN1.svg"), width = 2.8, height = 2)

#write numbers for Sarah
write_xlsx(neuts,here::here("04_Analysis", "data_objects", "paperfigures", "Figure 6", "XBBJN1Neut_ByBias.xlsx"))
######