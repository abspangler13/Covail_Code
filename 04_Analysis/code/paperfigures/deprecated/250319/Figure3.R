library(ggplot2)
library(dplyr)
library(here)
library(Seurat)
library(readxl)
library(iNEXT)
library(tidyseurat)

set.seed(1)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObj %>% filter(ClusterLabel != "Naive")

pIgG <- seuObj %>% tidyseurat::join_features(features = "P-IgG", assay = "Prot")
pIgG <- pIgG@meta.data %>% select(CELL, .feature, .abundance_Prot) %>% pivot_wider(names_from=.feature,values_from=.abundance_Prot)

df <- seuObj@meta.data %>% mutate(PIGG = pIgG$`P-IgG`[match(CELL, pIgG$CELL)],
                                  PIGG = PIGG + abs(min(PIGG)) + 1)

df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                    Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                    Booster == "Prototype" ~ "Prototype mRNA"))

#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#045275")

#####
#Include flow data
flow <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx")) %>% 
  filter(`Subject ID` %in% df$Subject & !duplicated(.))

#Prepare probe metrics
#Proto+Beta+
flow$`Proto+/Beta+/BA1+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+/BA1+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1 | Freq. of IgG` * 100
flow$`Proto+/Beta+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta | Freq. of IgG` * 100

#Proto Only
flow$`Proto+` <- flow$`Live/IgG/Proto+/Proto | Freq. of IgG` * 100
flow$`Proto+/BA1+/XBB+` <- flow$`Live/IgG/Proto+/Proto-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/BA1+` <- flow$`Live/IgG/Proto+/Proto-BA1 | Freq. of IgG` * 100
flow$`Proto+/XBB+` <- flow$`Live/IgG/Proto+/Proto-XBB | Freq. of IgG` * 100

#Beta Only
flow$`Beta+` <- flow$`Live/IgG/Beta+/Beta | Freq. of IgG` * 100
flow$`Beta+/BA1+/XBB+` <- flow$`Live/IgG/Beta+/Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Beta+/BA+` <- flow$`Live/IgG/Beta+/Beta-BA1 | Freq. of IgG` * 100
flow$`Beta+/XBB+` <- flow$`Live/IgG/Beta+/Beta-XBB | Freq. of IgG` * 100

#the m'crons
flow$`BA1+` <- flow$`Live/IgG/Proto neg Beta neg/BA1 | Freq. of IgG` * 100
flow$`XBB+` <- flow$`Live/IgG/Proto neg Beta neg/XBB | Freq. of IgG` * 100
flow$`BA1+/XBB+` <- flow$`Live/IgG/Proto neg Beta neg/BA1-XBB | Freq. of IgG` * 100

#Sum each of the totals:
#stage 1- we want a comparison between proto and omi-restricted compartments
flow$TotalRBD <- rowSums(flow[,c(55:69)])
flow$ProtoOmi <- rowSums(flow[,c(55, 56, 60, 61)])
#####

#####
#Figure 3B: Probe specificity by CITESeq
stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

# #write a file for Sarah
# write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "misc", "240328_SheetsForSarah","CITESeq_CrossReactivity_Calculations.csv"))

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+" & stats$Infection == "N",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_TotalCrossReactiveOnly.png"),width = 3.5, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#try a barplot
stats <- df %>% filter(Infection == "N") %>%
  mutate(adj.ProtoOmi = factor(adj.ProtoOmi, levels = c("Proto+Omi+", "Proto+Omi-", "Proto-Omi+"))) %>%
  group_by(Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(Subject, Timepoint, adj.ProtoOmi, fill = list(Proportion = 0, n = 0)) %>%
  mutate(OfficialBooster = df$OfficialBooster[match(Subject, df$Subject)]) %>%
  group_by(OfficialBooster, Timepoint, adj.ProtoOmi) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats[stats$OfficialBooster %in% c("Omicron BA.1 mRNA", "Prototype mRNA", "Prototype + Omicron BA.1 mRNA"),]) +
  geom_bar(aes(x=Timepoint, y=mean, fill = adj.ProtoOmi), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(15))+
  #scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  #scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#d7c49eff",
                               "Proto+Omi-" = "#343148ff",
                               "Proto-Omi+" = "#D64161FF"))+
  ylab("Proportion")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_ProtoOmicronCrossReactiveBarplot_CITESeq.png"),width = 3.3, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 3C: Correlation between proportion cross-reactive by flow and by CITESeq
flow$PropProtoOmi <- flow$ProtoOmi / flow$TotalRBD
flow$SubjTime <- paste(flow$`Subject ID`, flow$`Time point Guess`, sep="_")

stats <- df %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject),
         Timepoint = str_remove(Timepoint, "Day "),
         Timepoint = ifelse(Timepoint == "0", "1", Timepoint),
         SubjTime = paste(Subject, Timepoint, sep="_")) %>%
  filter(adj.ProtoOmi == "Proto+Omi+")

stats$PropFlow <- flow$PropProtoOmi[match(stats$SubjTime, flow$SubjTime)]

#plot the correlation
ggplot(stats, aes(x=Proportion, y = PropFlow))+
  geom_point(size = 2, shape = 21, aes(fill = OfficialBooster))+
  geom_abline(intercept=0, slope=1)+
  scale_fill_manual(values = allColors)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  ylab("Proportion Cross-Reactive By Flow")+
  xlab("Proportion Cross-Reactive By CITESeq")+
  ggtitle("Cross-Reactivity- CITESeq vs Flow")+
  theme_classic()+
  theme(plot.title = element_text(size=8), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 7),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqVsFlow_InfectedRemoved_PropCross.png"),width = 3.6, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()
#just per subject not timepoint
######

#####
#Figure 3d. Clonal overlap at all timepoints
nonSinglets <- unique(df$clone_subject_id[duplicated(df$clone_subject_id) | duplicated(df$clone_subject_id, fromLast=T)])
df$CloneStatus <- ifelse(df$clone_subject_id %in% nonSinglets, df$clone_subject_id, "Singlet")          

calcs <- df %>% filter(Infection == "N") %>%
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

stats <- df %>% filter(Infection == "N") %>%
          group_by(Subject, Timepoint, CloneStatusRefined) %>%
          summarize(
            n = n()) %>%
          mutate(Proportion = n / sum(n)) %>%
          ungroup() %>%
          complete(Subject, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
          mutate(OfficialBooster = df$OfficialBooster[match(Subject, df$Subject)]) %>%
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
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(15))+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3e_ClonalRelatednessOverTime.png"),width = 3.6, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()

#Try faceting by immunogen
stats <- df %>% filter(Infection == "N") %>%
  group_by(Subject, adj.ProtoOmi, Timepoint, CloneStatusRefined) %>%
  summarize(
    n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(Subject, adj.ProtoOmi, Timepoint, CloneStatusRefined, fill = list(Proportion = 0, n = 0)) %>%
  mutate(OfficialBooster = df$OfficialBooster[match(Subject, df$Subject)]) %>%
  group_by(OfficialBooster, adj.ProtoOmi, Timepoint, CloneStatusRefined) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = CloneStatusRefined), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(OfficialBooster), rows = vars(adj.ProtoOmi), labeller = label_wrap_gen(15))+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3e_ClonalRelatednessOverTime_byProbeSpecificity.png"),width = 3.7, height = 8.4, units = "in", device = "png", dpi = 600)
dev.off()

#####

#####
#Figure 3e: Simpson's index for clonal "evenness"
#let's try iNEXT
stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, adj.ProtoOmi, Timepoint, clone_id) %>%
  summarize(n = n()) %>%
  mutate(Chao1Richness = iNEXT::ChaoSimpson(n, datatype="abundance")$Estimator,
         Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))) %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, adj.ProtoOmi, Timepoint) %>%
  summarize(Chao1SimpsonDiversity = median(Chao1Richness))

#write a csv for sarah again
# write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "misc", "240328_SheetsForSarah","iNEXT_SpecificitySeparated_SimpsonDiversityEstimates.csv"))

#plot
ggplot(stats, aes(x = Timepoint, y = Chao1SimpsonDiversity))+
  #geom_errorbar(aes(ymin = median -sd, ymax = median + sd, color = OfficialBooster))+
  #geom_boxplot(aes(fill = OfficialBooster))+
  geom_line(aes(group = Subject, color = OfficialBooster), alpha = 0.5)+
  geom_point(aes(fill = OfficialBooster), shape=21, stroke = 0.3)+
  facet_grid(cols = vars(OfficialBooster), rows = vars(adj.ProtoOmi), labeller = label_wrap_gen(15))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  ylab("Estimated Simpson Index")+
  ylim(c(0.8, 1))+
  theme_classic()+
  theme(plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin=margin()),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_iNextClonalRichness_Chao1_ByProbeSpec.png"),width = 3.7, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#figure 3f: SHM overall
stats <- df %>% filter(Infection == "N") %>%
  group_by(OfficialBooster, Subject, Timepoint, mu_freq) %>%
  summarize(mean = mean(mu_freq)) %>%
  group_by(OfficialBooster, Timepoint) %>%
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_lineplot.png"),width = 2, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

#####

#####
#Figure 3g
#fold change
#plot fold change on group basis
stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+" & Infection == "N") %>%
                group_by(OfficialBooster, Subject, adj.ProtoOmi, Timepoint) %>%
                summarize(meanSHM = mean(mu_freq)) %>%
                mutate(FoldChange = meanSHM / meanSHM[1]) %>%
                group_by(OfficialBooster, Timepoint, adj.ProtoOmi) %>%
                summarize(meanFold = mean(FoldChange),
                          n= n(),
                          sd = sd(FoldChange)) %>%
                mutate(se = sd / sqrt(n))

ggplot(stats, aes(x = Timepoint, y= meanFold, fill = OfficialBooster))+
  geom_errorbar(aes(ymin = meanFold - se, ymax = meanFold + se, color = OfficialBooster), width = 0.4, alpha=0.9)+
  geom_line(aes(group = OfficialBooster, color = OfficialBooster))+
  geom_point(shape =21, size = 2, stroke = 0.3)+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Fold Change in SHM")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  geom_hline(yintercept = 1, linetype = "longdash")+
  facet_grid(cols = vars(adj.ProtoOmi)) +
  ylim(c(0.85, 1.3))+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMBySpec_FoldChange.png"),width = 2.5, height = 2, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 2H: lineages over time- supplemental
# calcs <- df %>% filter(Infection == "N") %>%
#           group_by(clone_subject_id, Timepoint) %>%
#           summarize(n = n()) %>%
#           mutate(AllTimes = length(intersect(Timepoint, c("Day 0", "Day 15", "Day 90", "Day 180"))) == 4)
# 
# stats <- df %>% filter(Infection == "N" & clone_subject_id %in% calcs$clone_subject_id[calcs$AllTimes == TRUE]) %>%
#           group_by(OfficialBooster, clone_subject_id, Timepoint) %>%
#           summarize(medSHM = median(mu_freq))
# 
# ggplot(stats, aes(x = Timepoint, y = medSHM))+
#   geom_line(aes(group = clone_subject_id, color = OfficialBooster), alpha= 0.5)+
#   geom_point(shape = 21, size =2, aes(fill = OfficialBooster))+
#   ggtitle("Median SHM")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(15))+
#   theme_classic()+
#   theme(title = element_text(size = 8),
#         axis.title.y = element_text(size=7),
#         axis.title.x = element_text(size=0),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=7),
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
#         panel.spacing = unit(0.3, "lines"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMByLineageOverTime.png"),width = 4, height = 2.2, units = "in", device = "png", dpi = 600)
# dev.off()
#####

#####
#Figure 3H alternative: violin plots for shm
check <- df %>% filter(Infection=="N")

ggplot(check, aes(x = Timepoint, y = mu_freq, fill = OfficialBooster))+
  # geom_point(shape =21, aes(fill = OfficialBooster), position = position_jitter(width = 0.2))+
  # stat_summary(geom="crossbar",
  #              fun = median,
  #              fun.min = median,
  #              fun.max = median,
  #              width = 0.8,
  #              linewidth = 0.4)+
  geom_violin()+
  ylab("Frequency of VH Mutations")+
  geom_boxplot(width=0.2, outlier.size = 0)+
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall.png"),width = 4, height = 2.2, units = "in", device = "png", dpi = 600)
dev.off()

#facet by immunogen
ggplot(check, aes(x = Timepoint, y = mu_freq, fill = OfficialBooster))+
  geom_violin()+
  ylab("Frequency of VH Mutations")+
  geom_boxplot(width=0.2, outlier.size = 0)+
  facet_grid(cols = vars(OfficialBooster), rows = vars(adj.ProtoOmi), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_facetedbyspec.png"),width = 4.2, height = 6.6, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 3: Clonal lineages over timee
lastingLins <- df %>% filter(Infection == "N") %>%
                group_by(clone_subject_id) %>%
                filter(length(unique(Timepoint)) == 4) %>%
                group_by(OfficialBooster, clone_subject_id, Timepoint) %>%
                summarize(medianSHM = median(mu_count))

ggplot(lastingLins, aes(x = Timepoint, y = medianSHM, fill = OfficialBooster))+
  geom_line(aes(group= clone_subject_id, color = OfficialBooster), alpha = 0.5)+
  geom_point(shape = 21, stroke = 0.5)+
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  ylab("Median VH Mutation Frequency")+
  #ggtitle("SHM Per Clonal Lineage")+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMPerClonalLineage.png"),width = 4, height = 2.2, units = "in", device = "png", dpi = 600)
dev.off()


#####
#Figure 3??: prototype and BA.1-binding strength per group
#calculate an "affinity" metric, normalized by IgG
# df$IgGNormProto <- df$`Proto-RBD-PE` / df$PIGG
# df$IgGNormBA1 <- df$`BA1-RBD-PE` / df$PIGG
# df$IgGNormXBB <- df$`XBB-RBD-no-fluor` / df$PIGG
# 
# #plot this
# ggplot(df[df$c_call %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4") & df$Infection == "Y",], aes(x = Timepoint, y = `BA1-RBD-PE`))+
#   geom_jitter(aes(fill = OfficialBooster), shape = 21)+
#   #geom_violin(aes(fill = OfficialBooster))+
#   facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   scale_y_log10()+
#   theme_classic()+
#   theme(title = element_text(size = 12),
#         axis.title.y = element_text(size=12),
#         axis.title.x = element_text(size=0),
#         axis.text.x = element_text(size=12,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=12),
#         legend.position = "none",
#         strip.background = element_blank(),
#         strip.text = element_text(size = 12, face = "bold", margin = margin()),
#         panel.spacing = unit(0.3, "lines"))







#####Code graveyard
#####
###Figure 3: Donut plots schism'ed by reactivity 
# splitDF <- split(df, df$adj.ProtoOmi)
# splitDF <- splitDF[2:3] #let's not do omicron-single positives
# for(j in 1:length(splitDF)){
#   
#   filename <- paste0(names(splitDF)[j],"_NussenzweigStyleDonuts.pdf")
#   tempDF <- splitDF[[j]]
#   
#   pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", filename))
#   for(i in unique(tempDF$Subject)[order(unique(tempDF$Subject))]){
#     placeholder <- tempDF[tempDF$Subject == i,]
#     grp <- unique(placeholder$Booster)
#     inf <- unique(placeholder$Infection)
#     
#     nonSinglets <- unique(placeholder$clone_subject_id[duplicated(placeholder$clone_subject_id) | duplicated(placeholder$clone_subject_id, fromLast=T)])
#     placeholder$CloneStatus <- ifelse(placeholder$clone_subject_id %in% nonSinglets, placeholder$clone_subject_id, "Singlet")          
#     
#     placeholder$Timepoint <- ifelse(placeholder$Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", placeholder$Timepoint)
#     
#     calcs <- placeholder %>%
#       group_by(CloneStatus, Timepoint) %>%
#       summarize(n = n()) %>%
#       mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
#                              #length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Day 0 Expanded",
#                              length(unique(Timepoint)) > 1 ~ "Expanded",
#                              TRUE ~ "Single Timepoint"))
#     
#     placeholder <- placeholder %>%
#       group_by(Timepoint, CloneStatus) %>%
#       summarize(n= n()) %>%
#       mutate(Proportion = n / sum(n),
#              Total = sum(n),
#              CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
#              adj.CloneStatus = case_when( CloneStatus == "Singlet" ~ "Singlet",
#                                           #CloneStatus %in% calcs$CloneStatus[calcs$lab == "Day 0 Expanded"] ~ "Day 0 Expanded",
#                                           CloneStatus %in% calcs$CloneStatus[calcs$lab == "Expanded"] ~ "Expanded",
#                                           CloneStatus %in% calcs$CloneStatus[calcs$lab == "Single Timepoint"] ~ "Single Timepoint"),
#              adj.CloneStatus = fct(adj.CloneStatus, levels = c("Day 0 Expanded", "Expanded", "Single Timepoint", "Singlet")),
#              Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90/180")))%>%
#       arrange(adj.CloneStatus)%>%
#       mutate(ymax = cumsum(Proportion),
#              ymin = c(0, head(ymax, n=-1)))
#     
#     #placeholder$CloneStatus <- factor(placeholder$CloneStatus, levels= c("Singlet","Expanded", "Day 0 Expanded", "Single Timepoint"))
#     label <- c()
#     label[1] <- unique(placeholder$Total[placeholder$Timepoint == "Day 0"])
#     label[2] <- unique(placeholder$Total[placeholder$Timepoint == "Day 15"])
#     label[3] <- unique(placeholder$Total[placeholder$Timepoint == "Day 90/180"])
#     time <- c("Day 0", "Day 15", "Day 90/180")
#     dat_text <- data.frame(label = label, Timepoint = time)
#     dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90/180"))
#     
#     
#     p <- ggplot(placeholder)+
#       geom_rect(color= "black", linewidth=0.2, aes(ymax=ymax, ymin=ymin, xmax=3.3, xmin=2.5, fill=adj.CloneStatus))+
#       coord_polar(theta="y")+
#       xlim(c(2,4))+
#       scale_fill_manual(values = c("Singlet" = "#FFFFFF",
#                                    #"Day 0 Expanded" = "#1d4f4b",
#                                    "Expanded" = "#40b0a7",
#                                    "Single Timepoint" = "#ffbe6a"))+
#       ggtitle(paste0("Group, Subject: ",i, " ", grp, ", Infection: ",inf, " ", names(splitDF)[j]))+
#       guides(fill = "none")+
#       facet_grid(cols=vars(Timepoint))+
#       theme_void()+
#       theme(plot.title = element_text(hjust=0.5),
#             panel.spacing = unit(-4.5, "lines"))+
#       geom_text(data = dat_text,
#                 mapping = aes(x=-Inf, y=-Inf, label = label),
#                 hjust = 0.5,
#                 vjust = 0.5)
#     
#     print(p)
#   }
#   dev.off()
#   rm(p)
# }
# rm(calcs)
#####

#####
#Figure 3D: Tabulated data using Alakazam/iNEXT
#there are a range of estimated diversity orders- q = {0,1,2}. As q increases, response evenness (whether one group proportionally expands) is weighted over richness (outright number of clonal groups)
#I think we don't care so much about richness- richness is more a factor of the PBMC sampling we have here- so I am tempted to use q=1 (Shannon index) or q=2 (Simpson Index)
#since we really care about whether or not a specific clonal group domineers a response, let's settle for q=2 right now (Simpson index)
# stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
#           group_by(OfficialBooster, Infection, InfectionRange, Subject, adj.ProtoOmi, Timepoint, clone_id) %>%
#           summarize(n = n()) %>%
#           mutate(index = alakazam::calcDiversity(n, q=2)) %>%
#           summarize(index = median(index)) %>%
#           ungroup() %>%
#           mutate(Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))
# 
# check <- stats %>%
#           group_by(OfficialBooster, Infection, InfectionRange, Timepoint, adj.ProtoOmi) %>%
#           summarise(median = median(index),
#                     sd = sd(index))
# 
# #show individual data (boxplots are getting really old)
# pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_AlakazamDiversity_SimpsonIndex_ByProbeSpec.pdf"), width = 8, height = 4)
# ggplot(stats[stats$Infection == "N",], aes(x = Timepoint, y = index, fill= adj.ProtoOmi))+
#   #geom_boxplot()+
#   geom_point(shape = 21)+
#   geom_line(aes(group = SubjectProto))+
#   facet_grid(cols = vars(OfficialBooster))+
#   scale_fill_manual(values = c("Proto+Omi+" = "#546a76",
#                                 "Proto+Omi-" = "#88a0a8"))+
#   ggtitle("Simpson's Index for Prototype-Restricted vs Cross-Reactive Cells")+
#   ylab("Simpson Index of Diversity")+
#   theme_classic()+
#   theme(legend.key.size = unit(0.4, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10),
#         strip.background = element_blank())
# dev.off()
# 
# #population median
# pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_AlakazamDiversity_SimpsonIndex_ByProbeSpec_Populationmedianandsd.pdf"), width = 8, height = 4)
# ggplot(check[check$Infection == "N",], aes(x = Timepoint, y = median, color= OfficialBooster))+
#   #geom_errorbar(aes(ymin = median -sd, ymax = median + sd, color = OfficialBooster))+
#   geom_point(aes(shape = OfficialBooster))+
#   geom_line(aes(group = OfficialBooster))+
#   facet_grid(cols = vars(adj.ProtoOmi))+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   ggtitle("Simpson's Index for Prototype-Restricted vs Cross-Reactive Cells")+
#   ylab("Simpson Index of Diversity")+
#   theme_classic()+
#   theme(legend.key.size = unit(0.4, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10),
#         strip.background = element_blank())
# dev.off()
# 
# #Diversity of infected group
# #show individual data (boxplots are getting really old)
# pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_AlakazamDiversity_SimpsonIndex_ByProbeSpec_Infected.pdf"), width = 8, height = 5)
# ggplot(stats[stats$Infection == "Y",], aes(x = Timepoint, y = index, fill= adj.ProtoOmi))+
#   geom_boxplot()+
#   #geom_point(shape = 21)+
#   #geom_line(aes(group = SubjectProto))+
#   facet_grid(cols = vars(OfficialBooster), rows = vars(InfectionRange))+
#   scale_fill_manual(values = c("Proto+Omi+" = "#546a76",
#                                "Proto+Omi-" = "#88a0a8"))+
#   ggtitle("Simpson's Index for Prototype-Restricted vs Cross-Reactive Cells - Infected")+
#   ylab("Simpson Index of Diversity")+
#   theme_classic()+
#   theme(legend.key.size = unit(0.4, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12, angle= 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10),
#         strip.background = element_blank())
# dev.off()
# 
# #population median
# pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_AlakazamDiversity_SimpsonIndex_ByProbeSpec_Populationmedianandsd_Infected.pdf"), width = 8, height = 5)
# ggplot(check[check$Infection == "Y",], aes(x = Timepoint, y = median, color= OfficialBooster))+
#   #geom_errorbar(aes(ymin = median -sd, ymax = median + sd, color = OfficialBooster))+
#   geom_point(aes(shape = OfficialBooster))+
#   geom_line(aes(group = OfficialBooster))+
#   facet_grid(cols = vars(adj.ProtoOmi), rows = vars(InfectionRange))+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   ggtitle("Simpson's Index for Prototype-Restricted vs Cross-Reactive Cells - Infected")+
#   ylab("Simpson Index of Diversity")+
#   theme_classic()+
#   theme(legend.key.size = unit(0.4, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10),
#         strip.background = element_blank())
# dev.off()

# #let's try iNEXT
# stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
#         group_by(OfficialBooster, Infection, InfectionRange, Subject, adj.ProtoOmi, Timepoint, clone_id) %>%
#         summarize(n = n()) %>%
#         mutate(Chao1Richness = iNEXT::ChaoSimpson(n, datatype="abundance")$Estimator,
#                Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))) %>%
#   group_by(OfficialBooster, Infection, InfectionRange, Subject, adj.ProtoOmi, Timepoint) %>%
#   summarize(Chao1SimpsonDiversity = median(Chao1Richness))
# 
# #write a csv for sarah again
# write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "misc", "240328_SheetsForSarah","iNEXT_SpecificitySeparated_SimpsonDiversityEstimates.csv"))
# 
# #plot
# ggplot(stats[stats$Infection == "N",], aes(x = Timepoint, y = Chao1Richness))+
#   #geom_errorbar(aes(ymin = median -sd, ymax = median + sd, color = OfficialBooster))+
#   #geom_boxplot(aes(fill = OfficialBooster))+
#   geom_line(aes(group = Subject, color = OfficialBooster), alpha = 0.5)+
#   geom_point(aes(fill = OfficialBooster), shape=21)+
#   facet_grid(cols = vars(OfficialBooster), rows = vars(adj.ProtoOmi))+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   ggtitle("Simpson's Index for Different MBC Spec.")+
#   ylab("Estimated Simpson Index")+
#   theme_classic()+
#   theme(legend.key.size = unit(0.4, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10),
#         strip.background = element_blank())
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_iNextClonalRichness_Chao1_ByProbeSpec.png"),width = 12, height = 6, units = "in", device = "png", dpi = 600)
# dev.off()
# #####

#####
#Make stats sheet for shm
#####
statistics <-  df %>% filter(Infection  == "N") %>%
  select(Booster, Timepoint, mu_freq)

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 3", "SHMOverTime.csv"))

#clonal groups instead
statistics <- df %>% filter(Infection == "N") %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(OfficialBooster, clone_subject_id, Timepoint) %>%
  summarize(medianSHM = median(mu_count)) %>%
  pivot_wider(names_from = Timepoint, values_from = medianSHM)

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "stats", "Figure 3", "SHMOverTime_ByLineages.csv"))
