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
#plot the denovo population
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
  facet_grid(cols = vars(Status), scale = "free_x", labeller = label_wrap_gen(10))+
  ylab("Number of mAbs")+
  xlab("Donor")+
  scale_fill_manual(values = c("FALSE" = "#C5A586", "TRUE" = "#EEE9E3"), labels = c("FALSE" = "Prototype-BA.1+", "TRUE" = "Prototype+BA.1+"))+
  theme_classic()+
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        legend.key.size = unit(0.8, "lines"),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots","paperfigures", "keystone", "PerDonor_DeNovoResponse.svg"), width = 6, height = 3)

####compare SHM
citeRat <- citeRat %>% filter(!is.na(IgGQuant))

ratpigSHM <- observedMutations(citeRat, sequenceColumn = "sequence_alignment",
                               germlineColumn = "germline_alignment",
                               regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)

ratpigSHM2 <- ratpigSHM %>% mutate(DeNovo = BindingPopulation %in% c("BA.1-specific", "P-B+X+J+", "P-B+X+J-", "P-B+X-J+", "P-B+X-J-")) %>%
  filter(!(CorrectedBoost == "Prototype" & `Infection Status` == "Uninfected" & DeNovo))

ggplot(ratpigSHM2, aes(x= DeNovo, y= mu_freq*100))+
  geom_violin(aes(fill= DeNovo))+
  geom_boxplot(aes(fill = DeNovo), width = 0.2)+
  ylab("% VH Mutation")+
  scale_fill_manual(values = c("TRUE" = "#C5A586", "FALSE" = "#EEE9E3"))+
  scale_x_discrete(limits = c("FALSE", "TRUE"), labels = c("Prototype+BA.1+", "Prototype-BA.1+"))+
  theme_classic()+
  theme(text = element_text(size = 14, color = "black"),
        axis.text.x = element_text(size = 16, angle = 45, hjust= 1, vjust = 1, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "Denovo_Vs_Recall_SHM.svg"), width = 3, height = 4.1)
#####

#####
#finally, plot binding vs neut
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
  geom_point(aes(fill = CorrectedBoost, group = `Infection Status`, color =CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 2, stroke = 1)+
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
  theme(text = element_text(size = 14, color = "black"),
        legend.key.size = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, color = "black"),
        legend.position = "none",
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "BindingComparison_byInfection.svg"), width = 5, height = 4.1)

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
  geom_point(aes(fill = CorrectedBoost, group = `Infection Status`, color = CorrectedBoost), position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.55), size = 2, stroke = 1)+
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
  guides(fill = guide_legend(override.aes = list(shape=21), title = "Boost"),
         shape = guide_legend(override.aes = list(fill = "black")))+
  theme(text = element_text(size = 14, color = "black"),
        legend.key.size = unit(0.3, "lines"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1, color = "black"),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "NeutsComparison_ByInfection.svg"),width = 5, height = 4.1)


#####