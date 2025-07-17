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
#Figure 3B: Probe specificity by CITESeq
stats <- df %>%
  group_by(OfficialBooster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) #let's include infected donors but only at days 0 and 15

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Prototype+/Omicron+")+
  ggtitle("Omicron Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_TotalCrossReactiveOnly.png"),width = 3.3, height = 2.3, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_TotalCrossReactiveOnly.svg"),width = 3.3, height = 2.3)
dev.off()

########Plot out proto-specific population
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ggtitle("Prototype-Specific")+
  ylab("Proportion Prototype+/Omicron-")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),,
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_ProtoSpec.png"),width = 3.3, height = 2.3, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_ProtoSpec.svg"),width = 3.3, height = 2.3)
dev.off()

#Ba.1 specific population
ggplot(stats[stats$adj.ProtoOmi == "Proto-Omi+",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ggtitle("Omicron-Specific")+
  ylab("Proportion Prototype-/Omicron+")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),,
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_OmicronSpec.png"),width = 3.3, height = 2.3, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqData_InfectedRemoved_OmicronSpec.svg"),width = 3.3, height = 2.3)
dev.off()

#write stats
#cross reactive
stats <- df %>%
  filter(!c(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(adj.ProtoOmi == "Proto+Omi+") %>%
  pivot_wider(id_cols = !c(n),names_from = Timepoint, values_from = Proportion) %>%
  select(OfficialBooster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`)
write_xlsx(stats, here::here("04_Analysis", 'data_objects', "paperfigures", "Figure 3", "CrossReactive_Prop_OverTime.xlsx"))

#proto spec
stats <- df %>%
  filter(!c(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(adj.ProtoOmi == "Proto+Omi-") %>%
  pivot_wider(id_cols = !c(n),names_from = Timepoint, values_from = Proportion) %>%
  select(OfficialBooster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`)
write_xlsx(stats, here::here("04_Analysis", 'data_objects', "paperfigures", "Figure 3", "PrototypeSpecific_Prop_OverTime.xlsx"))

#ba1 spec
stats <- df %>%
  filter(!c(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(adj.ProtoOmi == "Proto-Omi+") %>%
  pivot_wider(id_cols = !c(n),names_from = Timepoint, values_from = Proportion) %>%
  select(OfficialBooster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`)
write_xlsx(stats, here::here("04_Analysis", 'data_objects', "paperfigures", "Figure S3", "OmicronSpecific_Prop_OverTime.xlsx"))
#########

#####

#####
#Figure 3C: Correlation between proportion cross-reactive by flow and by CITESeq
flow$PropProtoOmi <- flow$ProtoOmi / flow$TotalRBD
flow$SubjTime <- paste(flow$`Subject ID`, flow$TimepointC, sep="_")

stats <- df %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject),
         SubjTime = paste(Subject, Timepoint, sep="_")) %>%
  filter(adj.ProtoOmi == "Proto+Omi+")

stats$PropFlow <- flow$PropProtoOmi[match(stats$SubjTime, flow$SubjTime)]

#plot the correlation
ggplot(stats, aes(x=Proportion, y = PropFlow))+
  geom_point(size = 1.8, shape = 21, aes(fill = Booster))+
  scale_fill_manual(values = immunogenColors)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  ylab("Proportion by Flow")+
  xlab("Proportion By CITESeq")+
  ggtitle("Proportion Cross-Reactive")+
  #geom_abline(linetype = 2, intercept = 0, slope = 1)+
  theme_classic()+
  theme(plot.title = element_text(size=10), 
        axis.title.y = element_text(size=9),
        axis.title.x = element_text(size=9),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        panel.spacing = unit(0.1, "lines"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_blank(),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqVsFlow_InfectedRemoved_PropCross.png"),width = 3.3, height = 2.1, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_CITESeqVsFlow_InfectedRemoved_PropCross.svg"),width = 3.3, height = 2.1, units = "in")
dev.off()

#calculate correlation
cor.test(stats$PropFlow, stats$Proportion, method = "pearson")
######

####################################
#from here-on, we will be looking exclusively at cross-reactive cells:
crossDF <- df %>% filter(adj.ProtoOmi == "Proto+Omi+", Infection == "N")

######
#Clonality barplots
#Figure 3d. Clonal overlap at all timepoints
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3e_ClonalRelatednessOverTime.png"),width = 4, height = 2.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3e_ClonalRelatednessOverTime.svg"),width = 4, height = 2.8, units = "in")
dev.off()
#######

######
#plot SHM over time
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHM_crossreactives.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHM_crossreactives.svg"),width = 3, height = 1.8, units = "in")
dev.off()

#do stats
shm %>% select(Timepoint, Booster, mu_freq) %>%
write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 3", "mu_freq_OverTime.xlsx"))
#####

# #####
# #Plot by fold change
# stats <- shm %>%
#   group_by(Booster, Subject, Timepoint) %>%
#   summarize(mean = mean(mu_freq)) %>%
#   group_by(Booster, Subject) %>% arrange(Timepoint) %>% mutate(mean = mean / mean[1]) %>%
#   group_by(Booster, Timepoint) %>%
#   summarize(mean2 = mean(mean),
#             sd = sd(mean),
#             n = n())%>%
#   mutate(se = sd / sqrt(n))
# 
# ggplot(stats, aes(x = Timepoint, y = mean2, fill= Booster))+
#   geom_errorbar(aes(ymin = mean2 - se, ymax = mean2 + se, color = Booster), width = 0.2, alpha=0.9)+
#   geom_line(aes(group = Booster, color = Booster), alpha = 0.8)+
#   geom_point(shape =21, size = 2, stroke = 0.3)+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   ggtitle("Fold Change in SHM")+
#   ylab("FC in Mutations")+
#   geom_hline(yintercept = 1, linetype = 2)+
#   ylim(c(0.8, 2))+
#   theme_classic()+
#   #guides(fill = guide_legend(override.aes = list(size=1.4)))
#   theme(title = element_text(size = 7),
#         axis.title.y = element_text(size=7),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=7),
#         legend.text = element_text(size = 6),
#         legend.title = element_blank())
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_lineplot.png"),width = 2.6, height = 1.6, units = "in", device = "png", dpi = 600)
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_lineplot.svg"),width = 2.6, height = 1.6, units = "in")
# dev.off()
# 
# #write a stats sheet
# stats <- shm %>%
#   group_by(Booster, Subject, Timepoint) %>%
#   summarize(mean = mean(mu_freq)) %>%
#   group_by(Booster, Subject) %>% 
#   arrange(Timepoint) %>% 
#   mutate(mean = mean / mean[1]) %>%
#   pivot_wider(names_from = Timepoint, values_from = mean) %>%
# write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 3", "FCInMeanSHM.xlsx"))
# #####

#####
#SHM by lineage- only show lineages that exist at all timepoints
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_by_clone.png"),width = 2.9, height = 1.6, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMOverall_by_clone.svg"),width = 2.9, height = 1.6)
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
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 3", "SHM_perClonalLineage.xlsx"))
#####

#######
#Clonal evolution
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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ClonalEvolution.png"), width = 1.6, height =3)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ClonalEvolution.svg"), width = 1.6, height =2.8)
dev.off()
######

######
#Time to go to Dunkin's because we are making some donuts
#BA.1: 4848544848
#BA.1 + Prototype: 5053564848
#Prototype: 4955534848
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
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ClonalDonuts.png"), width = 2, height = 2.5, unit = "in", dpi = 1500)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ClonalDonuts.svg"), width = 2, height = 2.5)

#####

######
#random sarah graph
stats <- df %>%
  group_by(adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         adj.ProtoOmi = factor(adj.ProtoOmi, levels = c("Proto+Omi+",
                                                        "Proto+Omi-","Proto-Omi+")))

ggplot(stats, aes(x = 1, y = Proportion))+
  geom_bar(stat = "identity", position = "stack", aes(fill = adj.ProtoOmi))+
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  ylab("Proportion")+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "TotalCITESeqBarplot_Specificities.png"), width = 2.1, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "TotalCITESeqBarplot_Specificities.svg"), width = 2.1, height = 2)
######

#####
######plot bulk shm among these lasting lineages
lastingLins <- crossDF %>% filter(Infection == "N") %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4)

lastingLins %>% filter(Timepoint %in% c("Day 90", "Day 180"), Booster == "BA.1 And Prototype") %>% mutate(mu_freq = mu_freq * 100) %>%
  ggplot(aes(x = Timepoint, y = mu_freq))+
  geom_boxplot(fill = "#1D75BC", outlier.shape = NA)+
  geom_point(shape =21, aes(fill = Booster), position = position_jitter(width = 0.1), fill = "#1D75BC")+
  ggtitle("Bivalent Boost, Lasting Lineages")+
  ylab("% VH Mutation")+
  scale_x_discrete(limits = c("Day 90", "Day 180"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        plot.title = element_text(hjust = 0.5))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Bivalent_LastingLins_Bulk.png"), width = 2.3, height = 2)

#do some testing
wilcox.test(lastingLins$mu_freq[lastingLins$Timepoint == "Day 90"],
            lastingLins$mu_freq[lastingLins$Timepoint == "Day 180"], paired = FALSE)
  
#########plot lasting lins, scaling alpha by presence in evolving lineages
lastingLins <- crossDF %>% filter(Infection == "N") %>% mutate(mu_freq = mu_freq * 100) %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(OfficialBooster, clone_subject_id, Timepoint) %>%
  summarize(medianSHM = median(mu_freq)) %>%
  mutate(clone_subject_id = paste0("s", clone_subject_id, "_1"),
         sig= case_when(clone_subject_id %in% evolving$clone_id ~ evolving$sig[match(clone_subject_id, evolving$clone_id)],
                        TRUE ~ FALSE),
         sig = ifelse(is.na(sig), FALSE, sig))

ggplot(lastingLins, aes(x = Timepoint, y = medianSHM, fill = OfficialBooster))+
  geom_line(aes(group= clone_subject_id, color = OfficialBooster, alpha = sig))+
  geom_point(shape = 21, stroke = 0.5, aes(alpha = sig))+
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  ylab("Mean % VH Mutation")+
  #ggtitle("SHM Per Clonal Lineage")+
  theme_classic()+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMPerClonalLineage_highlightedsig.png"),width = 4, height = 2.2, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMPerClonalLineage_highlightedsig.svg"),width = 4, height = 2.2, units = "in", device = "png", dpi = 600)
dev.off()

#######Test differences in mean using paired analysis from just day 90 to d180
lastingLins <- crossDF %>% filter(Infection == "N") %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(medianSHM = median(mu_freq)) %>% filter(Timepoint %in% c("Day 90", "Day 180"), Booster == "BA.1 And Prototype")

ggplot(lastingLins, aes(x = Timepoint, y = medianSHM*100, fill = Booster))+
  geom_line(aes(group= clone_subject_id, color = Booster), alpha = 0.5)+
  geom_point(shape = 21, stroke = 0.5)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 90", "Day 180"))+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  ylab("Median VH Mutation Frequency")+
  #ggtitle("SHM Per Clonal Lineage")+
  ylim(0,11)+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_ClonalLineageSHM_d90_d180_bivalent.png"),width = 2, height = 2.2, units = "in", device = "png", dpi = 600)
dev.off()

#test
lastingLins <- crossDF %>% filter(Infection == "N") %>%
  group_by(clone_subject_id) %>%
  filter(length(unique(Timepoint)) == 4) %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(medianSHM = median(mu_freq)) %>% filter(Timepoint %in% c("Day 90", "Day 180"), Booster == "BA.1 And Prototype") %>%
  pivot_wider(names_from = Timepoint, values_from = medianSHM)

wilcox.test(lastingLins$`Day 90`,lastingLins$`Day 180`, paired= TRUE)
t.test(lastingLins$`Day 90`, lastingLins$`Day 180`,paired= TRUE)

######

######
#do proportion of probe specificities over time
#Figure 3B: Probe specificity by CITESeq
stats <- df %>% filter(Infection == "N") %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) #let's include infected donors but only at days 0 and 15

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Prototype+/Omicron+")+
  ggtitle("Omicron Cross-Reactive")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "CITESeq_ProportionCrossReactive_NoInfected.svg"),width = 3.3, height = 2.3)
dev.off()

#stats
stats <- df %>%
  filter(Infection == "N") %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(adj.ProtoOmi == "Proto+Omi+") %>%
  pivot_wider(id_cols = !c(n),names_from = Timepoint, values_from = Proportion) %>%
  select(OfficialBooster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`)
write_xlsx(stats, here::here("04_Analysis", 'data_objects', "paperfigures", "Figure 3", "CrossReactive_Prop_OverTime_NoInfected.xlsx"))

########Prototype specifics
stats <- df %>% filter(Infection == "N") %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) #let's include infected donors but only at days 0 and 15


ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=Proportion, fill=OfficialBooster))+
  geom_line(alpha = 0.5, aes(group = Subject, color = OfficialBooster))+
  geom_point(shape=21, aes(fill=OfficialBooster), stroke = 0.3)+
  ylab("Proportion Prototype+/Omicron-")+
  ggtitle("Prototype+/Omicron-")+
  facet_grid(cols=vars(OfficialBooster), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic() +
  theme(plot.title = element_text(size = 8, hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "CITESeq_ProportionPrototype_NoInfected.svg"),width = 3.3, height = 2.3)
dev.off()

#stats
stats <- df %>%
  filter(Infection == "N") %>%
  group_by(OfficialBooster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(adj.ProtoOmi == "Proto+Omi-") %>%
  pivot_wider(id_cols = !c(n),names_from = Timepoint, values_from = Proportion) %>%
  select(OfficialBooster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`)
write_xlsx(stats, here::here("04_Analysis", 'data_objects', "paperfigures", "Figure 3", "PrototypeOnly_Prop_OverTime_NoInfected.xlsx"))


#####





















############suppl
#####
#Calculate diversity
# diversityDF <- df %>% filter(Infection == "N") %>%
#   filter(Timepoint %in% c("Day 0", "Day 15"), adj.ProtoOmi != "Proto-Omi+") %>%
#   mutate(BoosterDay = paste(OfficialBooster, Subject, Timepoint, adj.ProtoOmi))
# 
# ###method 1: calculating using the alphadiversity function
# curve <- alakazam::alphaDiversity(diversityDF, group="BoosterDay", min_q=2, max_q=2, min_n=20) #calculate bootstrapped diversity estimate at q=2 (Simpson)
# 
# calculatedDiversity <- curve@diversity %>% filter(q == 2) %>% #convert data to graphable form
#   mutate(Booster = str_remove(BoosterDay, " [0-9]+ Day [0-9]+ Proto[+-]Omi[+-]"),
#          Timepoint = str_extract(BoosterDay, "Day [0-9]+"),
#          adj.ProtoOmi = str_extract(BoosterDay, "Proto[+-]Omi[+-]"),
#          Donor = str_extract(BoosterDay, "(?<=mRNA )[0-9]+")
#          ) %>%
#   group_by(Donor, adj.ProtoOmi) %>%
#   mutate(AllTimes = length(unique(Timepoint)) < 2,
#          adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Omicron Cross-Reactive",
#                                   adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Specific")) %>%
#   filter(!AllTimes)
# 
# ggplot(calculatedDiversity, aes(x = Timepoint, y = d))+
#   geom_point(aes(fill = Booster), shape =21)+
#   geom_line(aes(group = Donor, color = Booster))+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+  
#   facet_grid(rows = vars(Booster), cols = vars(adj.ProtoOmi), axes = "all")+
#   ylim(8, 21)+
#   ylab("Simpson's Index")+
#   theme_classic()+
#   theme(text = element_text(size = 8),
#         strip.background = element_blank(),
#         strip.text = element_text(face = "bold"),
#         legend.position = "none")
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_diversity_diversityCurve_perdonor.png"),width = 3, height = 4, units = "in", device = "png", dpi = 1200)
# dev.off()
# 
# #write a table
# calculatedDiversity %>% select(!c(d_sd, d_lower, d_upper, e, e_lower, e_upper)) %>% select(Booster, Donor, Timepoint, adj.ProtoOmi, q, d) %>% pivot_wider(names_from =  Timepoint, values_from = d) %>%
# write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 3", "diversity_perindividual_bootstrap.xlsx"))
# 
# ####method 3: calculate curve in bulk?
# # diversityDF2 <- df %>%
# #   #filter(Timepoint %in% c("Day 0", "Day 15"), adj.ProtoOmi != "Proto-Omi+") %>%
# #   filter(adj.ProtoOmi != "Proto-Omi+") %>%
# #   mutate(BoosterDay = paste(OfficialBooster, Timepoint, adj.ProtoOmi))
# # curve <- alakazam::alphaDiversity(diversityDF2, group="BoosterDay", min_q=2, max_q=2, clone = "clone_subject_id") #calculate bootstrapped diversity estimate at q=2 (Simpson)
# # 
# # calculatedDiversity <- curve@diversity %>% filter(q == 2) %>% #convert data to graphable form
# #   mutate(Booster = str_remove(BoosterDay, " Day [0-9]+ Proto[+-]Omi[+-]"),
# #          Timepoint = str_extract(BoosterDay, "Day [0-9]+"),
# #          adj.ProtoOmi = str_extract(BoosterDay, "Proto[+-]Omi[+-]")
# #   )
# # 
# # ggplot(calculatedDiversity, aes(x = Timepoint, y = d))+
# #   geom_errorbar(aes(ymax = d_upper, ymin = d_lower, color = adj.ProtoOmi))+
# #   geom_line(aes(group = adj.ProtoOmi, color = adj.ProtoOmi))+
# #   geom_point(aes(fill = adj.ProtoOmi), shape =21)+
# #   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
# #   scale_fill_manual(values = c("Proto+Omi+" = "green",
# #                                "Proto+Omi-" = "orange"))+
# #   scale_color_manual(values = c("Proto+Omi+" = "green",
# #                                "Proto+Omi-" = "orange"))+
# #   ylab("Simpson's Diversity Index")+
# #   # scale_fill_manual(values = allColors)+
# #   # scale_color_manual(values = allColors)+  
# #   facet_grid(cols = vars(Booster))+
# #   theme_classic()+
# #   theme(text = element_text(size = 8),
# #         strip.background = element_blank(),
# #         axis.text.x = element_text(angle = 45, hjust = 1, vjust =1 ))
# # ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_diversity_diversityCurve.png"),width = 4, height = 2.5, units = "in", device = "png", dpi = 1200)
# # dev.off()
# 
# #just calculate average clone size over time??
# # diversityDF3 <- df %>% filter(Infection == "N") %>%
# #   filter(adj.ProtoOmi != "Proto-Omi+") %>%
# #   mutate(BoosterDay = paste(OfficialBooster, Subject, Timepoint, adj.ProtoOmi))
# # 
# # stats <- diversityDF3 %>%
# #         group_by(Booster, Subject, adj.ProtoOmi, Timepoint, clone_subject_id) %>%
# #         summarize( n = n()) %>%
# #         group_by(Booster, Subject, adj.ProtoOmi, Timepoint) %>%
# #         summarize(mean = mean(n)) %>%
# #         group_by(Booster, Subject, Timepoint) %>%
# #         mutate(meanRatio = mean / mean[adj.ProtoOmi == "Proto+Omi-"],
# #                adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Omicron Cross-Reactive",
# #                                         adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Specific"))
# # 
# # ggplot(stats, aes(x = Timepoint, y= mean))+
# #   geom_hline(yintercept = 1, linetype = 2)+
# #   geom_point(shape = 21, aes(fill = Booster))+
# #   geom_line(alpha = 0.5, aes(group = Subject, color = Booster))+
# #   facet_grid(rows = vars(Booster), cols = vars(adj.ProtoOmi), axes = "all")+
# #   ylab("Mean Clonal Lineage Size")+
# #   ylim(0.9,3.3)+
# #   scale_fill_manual(values = immunogenColors)+
# #   scale_color_manual(values = immunogenColors)+
# #   theme_classic()+
# #   theme(text = element_text(size = 8),
# #         strip.background = element_blank(),
# #         strip.text = element_text(face = "bold"),
# #         legend.position = "none")
# # ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_diversity_meanlineagessize.png"),width = 3, height = 4, units = "in", device = "png", dpi = 1200)
# # 
# # #write a table
# # stats %>% select(!meanRatio) %>% pivot_wider(id_cols = c(Booster, Subject, adj.ProtoOmi), values_from = mean, names_from = Timepoint) %>%
# # write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S3", "MeanClonalLineageSize.xlsx"))
# # #####
# 
# #####
# #Figure 3g
# #fold change
# #plot fold change on group basis
# stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+" & Infection == "N") %>%
#                 group_by(OfficialBooster, Subject, adj.ProtoOmi, Timepoint) %>%
#                 summarize(meanSHM = mean(mu_freq)) %>%
#                 mutate(FoldChange = meanSHM / meanSHM[1]) %>%
#                 group_by(OfficialBooster, Timepoint, adj.ProtoOmi) %>%
#                 summarize(meanFold = mean(FoldChange),
#                           n= n(),
#                           sd = sd(FoldChange)) %>%
#                 mutate(se = sd / sqrt(n))
# 
# ggplot(stats, aes(x = Timepoint, y= meanFold, fill = OfficialBooster))+
#   geom_errorbar(aes(ymin = meanFold - se, ymax = meanFold + se, color = OfficialBooster), width = 0.4, alpha=0.9)+
#   geom_line(aes(group = OfficialBooster, color = OfficialBooster))+
#   geom_point(shape =21, size = 2, stroke = 0.3)+
#   ylab("Fold Change")+
#   xlab("Timepoint")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   ggtitle("Fold Change in SHM")+
#   scale_fill_manual(values = allColors)+
#   scale_color_manual(values = allColors)+
#   geom_hline(yintercept = 1, linetype = "longdash")+
#   facet_grid(cols = vars(adj.ProtoOmi)) +
#   ylim(c(0.85, 1.3))+
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
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_SHMBySpec_FoldChange.png"),width = 2.5, height = 2, units = "in", device = "png", dpi = 600)
# dev.off()
# #####