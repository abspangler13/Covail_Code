library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(gridExtra)
library(viridisLite)
library(Seurat)
library(tidyseurat)
library(forcats)
library(ggprism)
library(readxl)
library(ggpubr)

#load in the seurat object
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObj %>% filter(ClusterLabel != "Naive")
df <- seuObj@meta.data

df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"))

#repeat same figures as figure 3 for uninfected group
#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#045275")

#specificity proportion
stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+" & stats$Infection == "Y",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), rows= vars(InfectionRange), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "Infected_CITESeqData_TotalCrossReactiveOnly.png"),width = 3.5, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#let's try just omicron-binding broadly (i.e. cross-reactivity + omicron-specificity)
df$OmicronBinding <- ifelse(df$adj.ProtoOmi %in% c("Proto+Omi+", "Proto-Omi+"), "Omicron+", "Omicron-")

stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, OmicronBinding) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

ggplot(stats[stats$OmicronBinding == "Omicron+" & stats$Infection == "Y",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), rows= vars(InfectionRange), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "Infected_CITESeqData_TotalCrossReactiveOnly.png"),width = 3.5, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()


#Clonal overlap
nonSinglets <- unique(df$clone_subject_id[duplicated(df$clone_subject_id) | duplicated(df$clone_subject_id, fromLast=T)])
df$CloneStatus <- ifelse(df$clone_subject_id %in% nonSinglets, df$clone_subject_id, "Singlet")          

calcs <- df %>% filter(Infection == "Y") %>%
  group_by(CloneStatus, Timepoint) %>%
  summarize(n = n()) %>%
  mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                         length(unique(Timepoint)) > 1 ~ "Expanded",
                         TRUE ~ "Single Timepoint"))

df$CloneStatusRefined <- calcs$lab[match(df$CloneStatus, calcs$CloneStatus)]
df$CloneStatusRefined <- factor(df$CloneStatusRefined, levels = c("Expanded Pre-Vax",
                                                                  "Expanded",
                                                                  "Single Timepoint",
                                                                  "Singlet"))

o <- unique(df$Subject[df$Booster == "Omicron"])
p <- unique(df$Subject[df$Booster == "Prototype"])
op <- unique(df$Subject[df$Booster == "Omicron And Prototype"])
fifNine <- unique(df$Subject[df$InfectionRange == "Between Days 15-90"])
nineOne <- unique(df$Subject[df$InfectionRange == "Between Days 90-180"])

stats <- df %>% filter(Infection == "Y") %>%
  group_by(InfectionRange, Booster, Subject, Timepoint, CloneStatusRefined) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(InfectionRange, Booster, Subject, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
  filter(
    (Subject %in% o & Booster == "Omicron") | (Subject %in% p & Booster == "Prototype") | (Subject %in% op & Booster == "Omicron And Prototype"),
    (Subject %in% fifNine & InfectionRange == "Between Days 15-90") | (Subject %in% nineOne & InfectionRange == "Between Days 90-180")
  ) %>%
  mutate(OfficialBooster = df$OfficialBooster[match(Subject, df$Subject)]) %>%
  group_by(InfectionRange, Booster, Timepoint, CloneStatusRefined) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = CloneStatusRefined), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(Booster), rows = vars(InfectionRange), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Expanded Pre-Vax" = "#4558A7",
                               "Expanded" = "#7598c0",
                               "Single Timepoint" = "#EBBB59",
                               "Singlet" = "#F2F3F4"))+
  ylab("Proportion")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_text(size = 0),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "InfectedClonalOverlap.png"),width = 3.6, height = 5.6, units = "in", device = "png", dpi = 600)
dev.off()

#lineplot of mean shm
stats <- df %>% filter(Infection == "Y") %>%
  group_by(OfficialBooster, InfectionRange, Subject, Timepoint, mu_freq) %>%
  summarize(mean = mean(mu_freq)) %>%
  group_by(OfficialBooster, InfectionRange, Timepoint) %>%
  summarize(mean2 = mean(mean),
            sd = sd(mean),
            n = n())%>%
  mutate(se = sd / sqrt(n))

ggplot(stats, aes(x = Timepoint, y = mean2, fill= OfficialBooster))+
  geom_errorbar(aes(ymin = mean2 - se, ymax = mean2 + se, color = OfficialBooster), width = 0.4, alpha=0.9)+
  geom_line(aes(group = OfficialBooster, color = OfficialBooster), alpha = 0.8)+
  geom_point(shape =21, size = 2, stroke = 0.3)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(cols = vars(InfectionRange))+
  ggtitle("Mean SHM Per Group")+
  ylab("Frequency of VH Mutations")+
  ylim(c(0,0.075))+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "Infected_SHMOverall_lineplot.png"),width = 4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

####let's look at SHM over time on a clonal lineage basis
stats <- df %>%
          filter(Infection == "Y") %>%
          group_by(OfficialBooster, InfectionRange, clone_subject_id) %>%
          filter(!(length(intersect(Timepoint, c("Day 0", "Day 15", "Day 90", "Day 180"))) < 4)) %>%
          group_by(OfficialBooster, InfectionRange, clone_subject_id, Timepoint) %>%
          summarize(mean = mean(mu_freq))

ggplot(stats, aes(x = Timepoint, y = mean, fill = OfficialBooster))+
  geom_line(aes(group = clone_subject_id, color = OfficialBooster), alpha=0.5)+
  geom_point(shape =21, size = 2, stroke = 0.2)+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols = vars(OfficialBooster), rows = vars(InfectionRange))+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "infected_clonalineages_overtime.png"),width = 4.5, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()

####let's compare infected and uninfected more directly
stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion))+
  geom_line(alpha = 0.5, aes(group = Subject, color = InfectionRange))+
  geom_point(shape=21, aes(fill=InfectionRange), stroke = 0.3)+
  ylab("Proportion Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  #scale_fill_manual(values = allColors)+
  #scale_color_manual(values = allColors)+
  theme_classic() +
  theme(axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        panel.spacing = unit(0.1, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "Infected_CITESeqData_TotalCrossReactiveOnly.png"),width = 3.5, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#####let's look again
#####
#Let's try a population bar chart over time
#make a pivoted table
stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  group_by(OfficialBooster, InfectionRange, Timepoint, adj.ProtoOmi) %>%
  summarize(mean = mean(Proportion),
            sd = sd(Proportion),
            n= n()) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum)) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), 0, .)))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = adj.ProtoOmi), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=1-cumusum-se, ymax= ifelse(1-cumusum+se >= 1, 1, 1-cumusum+se)), width=0.2, alpha=1) +
  facet_grid(cols = vars(OfficialBooster), rows = vars(InfectionRange), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  # scale_fill_manual(values = c("Proto-BA.1-" = "#FFB000",
  #                              "Proto-BA.1+" = "#FE6100",
  #                              "Proto+BA.1-" = "#DC267F",
  #                             "Proto+BA.1+" = "#648FFF"))+
  scale_fill_brewer(palette = "Spectral")+
  ylab("Proportion")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_text(size = 0),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S1", "Figures1_ProbePositivePopulationsOverTime_ProtoOmi.png"),width = 3.3, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()
#####

####look at SHM by per donor basis
stats <- df %>%
  group_by(OfficialBooster, InfectionRange, Subject, Timepoint) %>%
  summarize(mean = mean(mu_freq),
            sd = sd(mu_freq),
            n = n())%>%
  mutate(se = sd / sqrt(n))

ggplot(stats, aes(x = Timepoint, y = mean, fill= OfficialBooster))+
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, color = OfficialBooster), width = 0.4, alpha=0.9)+
  geom_line(aes(group = Subject, color = OfficialBooster), alpha = 0.8)+
  geom_point(shape =21, size = 2, stroke = 0.3)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(cols = vars(OfficialBooster),rows = vars(InfectionRange))+
  ggtitle("Mean SHM Per Group")+
  ylab("Frequency of VH Mutations")+
  ylim(c(0,0.075))+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "infectedfigures", "Infected_SHMOverall_lineplot.png"),width = 4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()



#check IgA
stats <- df %>%
          group_by(OfficialBooster, InfectionRange, Timepoint, c_call) %>%
          summarize(n = n()) %>%
          mutate(Prop = n / sum(n))

ggplot(stats, aes(fill=c_call, y=Prop, x=Timepoint))+
  geom_bar(position="stack",stat="identity", color="black")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_brewer(palette = "Spectral")+
  facet_grid(rows = vars(InfectionRange), cols=(vars(OfficialBooster)))+
  ylab("Proportion of all RBD+ Cells")+
  xlab("Timepoint")+
  theme_classic()+
  theme(axis.text.x=element_text(size=12, angle = 90), axis.text.y=element_text(size=12))



###let's plot the distributions of SHM for each group
infected <- df %>% filter(Infection == "Y") %>% mutate(Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))
infected2 <- df %>% mutate(Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))

ggplot(infected2, aes(x = mu_count, fill = InfectionRange))+
  geom_histogram()+
  facet_grid(cols = vars(Timepoint), rows = vars(Infection))+
  geom_vline(xintercept = 10)+
  theme_classic()+
  theme()




##tabulated SHM
stats <- df %>%
          filter(ClusterLabel != "Naive") %>%
          group_by(InfectionRange, Timepoint) %>%
          summarize(medianSHM = median(mu_count),
                    meanSHM = mean(mu_count))
