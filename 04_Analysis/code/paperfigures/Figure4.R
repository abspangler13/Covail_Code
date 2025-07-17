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
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")))

#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#2AB673", 
               "Prototype + Omicron BA.1 mRNA" = "#1D75BC",
               "Prototype mRNA" = "#FBB042")

immunogenColors <- c("Prototype" = "#FBB042",
                     "BA.1 And Prototype" = "#1D75BC",
                     "BA.1" = "#2AB673")
#####read in flow data
flowRaw <- read_xlsx(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Unfiltered_COVAILFlowDataset_250318.xlsx"))

flow <- flowRaw %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         Immunogen = str_replace_all(Immunogen, " \\+ ", "/"),
         Booster = str_replace_all(Booster, " \\+ ", "/"),
         TimepointC = case_when(Timepoint == 1 ~ "Day 0",
                                Timepoint == 15 ~ "Day 15",
                                Timepoint == 90 ~ "Day 90",
                                Timepoint == 180 ~ "Day 180"),
         TimepointC = factor(TimepointC, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))) %>% filter(Dataset != "Delta Panel")

#read in evolution data from Abby
evolving <- read.csv(here::here("01_raw-data", "Evo_dat_Timepoint_uniform.csv"))
#####

#####
#Make donut plots of clonality for each donor
crossDF <- df %>% filter(adj.ProtoOmi == "Proto+Omi+", Infection == "N")

crossDF2 <- crossDF %>% mutate(Timepoint = case_when(Timepoint %in% c("Day 90", "Day 180") ~ "Days 90/180",
                                                     TRUE ~ Timepoint)) %>%
  filter(Subject %in% c("4848544848", "5053564848", "4955534848"))

plotList <- list()
for(i in unique(crossDF2$Subject)){
  yuh <- crossDF2 %>% filter(Subject == i)
  
  nonSinglets <- unique(yuh$clone_subject_id[duplicated(yuh$clone_subject_id) | duplicated(yuh$clone_subject_id, fromLast=T)])
  yuh$CloneStatus <- ifelse(yuh$clone_subject_id %in% nonSinglets, yuh$clone_subject_id, "Singlet")
  
  calcs <- yuh %>% filter(Infection == "N") %>%
    group_by(CloneStatus, Timepoint) %>%
    summarize(n = n()) %>%
    mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                           length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                           length(unique(Timepoint)) > 1 ~ "Expanded",
                           TRUE ~ "Single Timepoint"))
  
  yuh$CloneStatusRefined <- calcs$lab[match(yuh$CloneStatus, calcs$CloneStatus)]
  yuh$CloneStatusRefined <- factor(yuh$CloneStatusRefined, levels = c("Expanded Pre-Vax",
                                                                      "Expanded",
                                                                      "Single Timepoint",
                                                                      "Singlet"))
  
  calcs <- yuh %>%
    group_by(Booster, Subject, Timepoint, CloneStatus) %>%
    summarize(n= n()) %>%
    mutate(Proportion = n / sum(n),
           Total = sum(n),
           ClonalOverlap = yuh$CloneStatusRefined[match(CloneStatus, yuh$CloneStatus)],
           #ClonalOverlap = ifelse(is.na(ClonalOverlap), "Singlet", ClonalOverlap),
           ClonalOverlap = factor(ClonalOverlap, levels = c("Expanded Pre-Vax",
                                                            "Expanded",
                                                            "Single Timepoint",
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
  
  p <- ggplot(calcs)+
    geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=ClonalOverlap))+
    coord_polar(theta="y")+
    xlim(c(2,4))+
    scale_fill_manual(values = c("Expanded Pre-Vax" = "#4558A7",
                                 "Expanded" = "#7598c0",
                                 "Single Timepoint" = "#EBBB59",
                                 "Singlet" = "#F2F3F4"))+
    ggtitle(unique(calcs$Booster))+
    guides(fill = "none")+
    geom_text(data = dat_text,
              mapping = aes(x=-Inf, y=-Inf, label = label),
              hjust = 0.5,
              vjust = 0.5,
              size = 2)+
    facet_grid(cols=vars(Timepoint), labeller = label_wrap_gen(1))+
    theme_void()+
    theme(plot.title = element_text(size = 3, hjust = 0.5, vjust = 0.5),
          strip.text.x = element_text(size = 6),
          panel.spacing = unit(-0.49999, "lines"))
  plotList[[i]] <- p
}

plotList <- plotList[c("4955534848", "5053564848", "4848544848")]
g <- arrangeGrob(grobs = plotList, nrow = length(plotList), ncol=1)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ClonalDonuts.png"), width = 2, height = 2.5, unit = "in", dpi = 1500)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ClonalDonuts.svg"), width = 2, height = 2.5)
#####

#####
#Make clonal barplots over time for cross-reactive cells
crossDF <- df %>% filter(adj.ProtoOmi == "Proto+Omi+", Infection == "N")

nonSinglets <- unique(crossDF$clone_subject_id[duplicated(crossDF$clone_subject_id) | duplicated(crossDF$clone_subject_id, fromLast=T)])
crossDF$CloneStatus <- ifelse(crossDF$clone_subject_id %in% nonSinglets, crossDF$clone_subject_id, "Singlet")

calcs <- crossDF %>% filter(Infection == "N") %>%
  group_by(CloneStatus, Timepoint) %>%
  summarize(n = n()) %>%
  mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                         length(unique(Timepoint)) > 1 ~ "Expanded",
                         TRUE ~ "Single Timepoint"))

crossDF$CloneStatusRefined <- calcs$lab[match(crossDF$CloneStatus, calcs$CloneStatus)]
crossDF$CloneStatusRefined <- factor(crossDF$CloneStatusRefined, levels = c("Expanded Pre-Vax",
                                                                            "Expanded",
                                                                            "Single Timepoint",
                                                                            "Singlet"))

stats <- crossDF %>% filter(Infection == "N") %>%
  group_by(Subject, Timepoint, CloneStatusRefined) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(Subject, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
  mutate(OfficialBooster = crossDF$OfficialBooster[match(Subject, crossDF$Subject)]) %>%
  group_by(OfficialBooster, Timepoint, CloneStatusRefined) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = CloneStatusRefined), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(12))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Expanded Pre-Vax" = "#4558A7",
                               "Expanded" = "#7598c0",
                               "Single Timepoint" = "#EBBB59",
                               "Singlet" = "#F2F3F4"))+
  ylab("Proportion")+
  theme_classic()+
  theme(legend.key.size = unit(0.3, 'cm'),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold", margin = margin()),
        panel.spacing = unit(0.35, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4e_ClonalRelatednessOverTime.png"),width = 4, height = 2.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4e_ClonalRelatednessOverTime.svg"),width = 4, height = 2.8, units = "in")
dev.off()
#####

#####
#make clone barplots for proto-specific cells
protoDF <- df %>% filter(adj.ProtoOmi == "Proto+Omi-", Infection == "N")

nonSinglets <- unique(protoDF$clone_subject_id[duplicated(protoDF$clone_subject_id) | duplicated(protoDF$clone_subject_id, fromLast=T)])
protoDF$CloneStatus <- ifelse(protoDF$clone_subject_id %in% nonSinglets, protoDF$clone_subject_id, "Singlet")

calcs <- protoDF %>% filter(Infection == "N") %>%
  group_by(CloneStatus, Timepoint) %>%
  summarize(n = n()) %>%
  mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                         length(unique(Timepoint)) > 1 ~ "Expanded",
                         TRUE ~ "Single Timepoint"))

protoDF$CloneStatusRefined <- calcs$lab[match(protoDF$CloneStatus, calcs$CloneStatus)]
protoDF$CloneStatusRefined <- factor(protoDF$CloneStatusRefined, levels = c("Expanded Pre-Vax",
                                                                            "Expanded",
                                                                            "Single Timepoint",
                                                                            "Singlet"))

stats <- protoDF %>% filter(Infection == "N") %>%
  group_by(Subject, Timepoint, CloneStatusRefined) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(Subject, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
  mutate(OfficialBooster = crossDF$OfficialBooster[match(Subject, crossDF$Subject)]) %>%
  group_by(OfficialBooster, Timepoint, CloneStatusRefined) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = CloneStatusRefined), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(12))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Expanded Pre-Vax" = "#4558A7",
                               "Expanded" = "#7598c0",
                               "Single Timepoint" = "#EBBB59",
                               "Singlet" = "#F2F3F4"))+
  ylab("Proportion")+
  theme_classic()+
  theme(legend.key.size = unit(0.3, 'cm'),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold", margin = margin()),
        panel.spacing = unit(0.35, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 8),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4e_ClonalRelatednessOverTime.png"),width = 4, height = 2.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4e_ClonalRelatednessOverTime.svg"),width = 4, height = 2.8, units = "in")
dev.off()
#####

#####
#SHM over time
shm <- crossDF %>% mutate(mu_freq = mu_freq * 100)

ggplot(shm, aes(x = Timepoint, y = mu_freq, fill = OfficialBooster))+
  geom_violin()+
  ylab("% VH Mutation")+
  geom_boxplot(width=0.2, outlier.size = 0)+
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  theme_classic()+
  ylim(0, 12.5)+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_SHM_crossreactives.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_SHM_crossreactives.svg"),width = 3, height = 1.8, units = "in")
dev.off()

#do stats
shm %>% select(Timepoint, Booster, mu_freq) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 4", "mu_freq_OverTime.xlsx"))
#####

#####
#Look on a per-lineage basis
clones <- shm %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  group_by(clone_subject_id) %>%
  mutate(Present = length(unique(Timepoint)) == 4) %>%
  filter(Present)

ggplot(clones, aes(x = Timepoint, y = mean, fill = Booster))+
  geom_line(aes(color = Booster, group = clone_subject_id), alpha= 0.5, linewidth = 0.8)+
  geom_point(shape = 21)+
  ylab("Mean % VH Mutation")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols= vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 7),
        strip.text = element_text(size = 8),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_SHMOverall_by_clone.png"),width = 2.9, height = 1.6, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_SHMOverall_by_clone.svg"),width = 2.9, height = 1.6)
dev.off()

#write a stats sheet
stats <- shm %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  group_by(clone_subject_id) %>%
  mutate(Present = length(unique(Timepoint)) == 4) %>%
  filter(Present) %>%
  group_by(Booster,clone_subject_id)%>% select(!c(n, Present)) %>%
  pivot_wider(names_from = Timepoint, values_from = mean) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 4", "SHM_perClonalLineage.xlsx"))
#####

#####
#make Abby's evolution plot
evolvingUninfected <- evolving %>% filter(!is.na(sig), Infection == "N", adj.ProtoOmi != "Proto-Omi+") %>% mutate(Booster = str_replace(Booster, "Omicron", "BA.1"))

summary <- evolvingUninfected %>%
  group_by(Booster) %>%
  summarize(n = n())

ggplot(evolvingUninfected)+
  geom_jitter(aes(fill = adj.ProtoOmi, alpha = sig, shape = adj.ProtoOmi, x = Booster, y= slope), width = 0.2)+
  scale_shape_manual(values = c("Proto+Omi+" = 21, "Proto+Omi-" = 22, "Proto-Omi+"= 24))+
  scale_fill_manual(values = c("Proto+Omi+" = "#FFA630", "Proto+Omi-" =  "#4DA1A9", "Proto-Omi+" = "#D7E8BA"))+
  #scale_y_continuous(position = "right")+
  ylab("Slope")+
  scale_alpha_discrete(guide = "none")+
  scale_x_discrete(limits= c("Prototype", "BA.1 And Prototype", "BA.1"))+
  ggtitle("Lineage Evolution")+
  geom_text(data = summary, mapping = aes(label = n, x = Booster, y = 0.0005), size = 3.5)+
  theme_classic()+
  guides(shape = guide_legend(nrow = 2))+
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ClonalEvolution.png"), width = 1.6, height =3)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ClonalEvolution.svg"), width = 1.6, height =2.8)
dev.off()
#####
