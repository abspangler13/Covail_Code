library(Seurat)
library(dplyr)
library(tidyseurat)
library(ggplot2)
library(ggalluvial)
library(readxl)
library(stringr)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))
seuObj <- seuObj %>% filter(ClusterLabel != "Naive")
df <- seuObj@meta.data
df$Timepoint <- factor(df$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

#####
#are resting mmemory and activated different as far as probe binding goes?
df$ActiveStatus <- case_when(df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)") ~ "Activated",
                             df$ClusterLabel %in% c("Resting IgG Memory", "Resting IgG Memory 2") ~ "Resting",
                             TRUE ~ "Irrelevant")

ggplot(df[df$Booster == "Omicron" & df$ActiveStatus != "Irrelevant" & df$Timepoint == "Day 15",], aes(x= Proto.RBD.PE, y=BA1.RBD.PE))+
  geom_point(shape =21, aes(fill = ProtoOmi))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  facet_grid(cols = vars(ActiveStatus))+
  ggtitle("Omicron Booster, AM1/2 vs Resting IgG Memories")+
  theme_classic()+
  theme()
#####

#####Trying the first indeterminate zone (BA.1 = 2-8)
#let's visualize the BA.1 plot and set a cutoff signal that we want for uncertain calls- I believe this is the 
ggplot(df[df$Booster == "Omicron" & df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"),], aes(x= Proto.RBD.PE, y=BA1.RBD.PE))+
  geom_point(shape =21, aes(fill = ProtoOmi))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  ggtitle("AM1/AM2 Cells for Omicron Booster")+
  facet_grid(cols = vars(Timepoint))+
  theme_classic()+
  geom_hline(yintercept = 8)+
  geom_hline(yintercept = 2)+
  theme() #it looks like from BA.1 signal 3-8, there are a decent number of cross-reactive cells intermingled with proto-specific

#for now, we'll label anything with a BA.1 signal of 3-8 as indeterminate to see if
#this has an outsized effect on BA.1 binders
df$AdjustedLabel <- case_when(df$BA1.RBD.PE >= 2 & df$BA1.RBD.PE <=8 ~ "Indeterminate",
                              TRUE ~ df$adj.ProtoOmi)

df$BA1.RBD.PE_Positive <- case_when(df$BA1.RBD.PE >= 2 ~ TRUE,
                              TRUE ~ df$BA1.RBD.PE_Positive)

df$ProtoOmi <- case_when(df$Proto.RBD.PE_Positive & (df$BA1.RBD.PE_Positive | df$XBB.RBD.no.fluor_Positive) ~ "Proto+Omi+",
                                !df$Proto.RBD.PE_Positive & (df$BA1.RBD.PE_Positive | df$XBB.RBD.no.fluor_Positive) ~ "Proto-Omi+",
                                df$Proto.RBD.PE_Positive & !(df$BA1.RBD.PE_Positive | df$XBB.RBD.no.fluor_Positive) ~ "Proto+Omi-",
                                TRUE ~ "Proto-Omi-")

adj.spec.P <- as.data.frame(df %>%
                              group_by(clone_subject_id,ProtoOmi) %>% 
                              dplyr::summarise(Freq = n()) %>% 
                              pivot_wider(names_from = ProtoOmi, values_from = Freq) %>% 
                              replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$clone_subject_id
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P$adj.ProtoOmi<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("adj.ProtoOmi","clone_subject_id")]
df <- df %>% left_join(adj.spec.P,by="clone_subject_id")

#let's see what group gets impacted most by this label (as far as activated memory goes)
ggplot(df[df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)") & df$Timepoint == "Day 15",], aes(x=Booster, y=1, fill = AdjustedLabel))+
  geom_col(position = "fill")+
  ggtitle("BA.1 Signal = 2-8 as indeterminate - AM1/2 Cells At Day 15")+
  theme_classic()+
  theme()

ggplot(df[df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)") & df$Timepoint == "Day 15",], aes(x=Booster, y=1, fill = adj.ProtoOmi.x))+
  geom_col(position = "fill")+
  ggtitle("Am1/Am2 Cells At Day 15 Prior to Dropping The Threshold")+
  theme_classic()+
  theme()

ggplot(df[df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)") & df$Timepoint == "Day 15",], aes(x=Booster, y=1, fill = adj.ProtoOmi.y))+
  geom_col(position = "fill")+
  ggtitle("BA.1 Signal = 2-8 as Positive - AM1/2 Cells At Day 15 With Clone Correction")+
  theme_classic()+
  theme()
#####

#####
#What prop of indet cells were cross-reactive?
ggplot(df[df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)") & df$Timepoint == "Day 15" & df$AdjustedLabel == "Indeterminate",], aes(x=Booster, y=1, fill = adj.ProtoOmi.x))+
  geom_col(position = "fill")+
  ggtitle("BA.1 Signal 2-8 set as indeterminate")+
  theme_classic()+
  theme()
#####

#####
#How does this change how this agrees w flow data?
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

#do a correlation
flow$PropProtoOmi <- flow$ProtoOmi / flow$TotalRBD
flow$SubjTime <- paste(flow$`Subject ID`, flow$`Time point Guess`, sep="_")

stats <- df %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi.y) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject),
         Timepoint = str_remove(Timepoint, "Day "),
         Timepoint = ifelse(Timepoint == "0", "1", Timepoint),
         SubjTime = paste(Subject, Timepoint, sep="_")) %>%
  filter(adj.ProtoOmi.y == "Proto+Omi+")

stats$PropFlow <- flow$PropProtoOmi[match(stats$SubjTime, flow$SubjTime)]

#plot the correlation
allColors <- c("Omicron" = "#7C1D6f", 
               "Omicron And Prototype" = "#DC3977",
               "Prototype" = "#045275")

ggplot(stats, aes(x=Proportion, y = PropFlow))+
  geom_point(size = 2, shape = 21, aes(fill = Booster))+
  geom_abline(intercept=0, slope=1)+
  scale_fill_manual(values = allColors)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  ylab("Proportion Cross-Reactive By Flow")+
  xlab("Proportion Cross-Reactive By CITESeq")+
  ggtitle("Cross-Reactive Response Comparison - 2-8")+
  theme_classic()+
  theme(plot.title = element_text(size=8), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "06_repertoire_analysis", "tinkering", "CITESeq_Vs_Flow_SettingIndeterminateRegionAsCrossReactive.png"),width = 3, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#stats thingy
stats <- df %>% filter(AdjustedLabel != "Indeterminate") %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi.x) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject),
         Timepoint = str_remove(Timepoint, "Day "),
         Timepoint = ifelse(Timepoint == "0", "1", Timepoint),
         SubjTime = paste(Subject, Timepoint, sep="_")) %>%
  filter(adj.ProtoOmi.x == "Proto+Omi+")

stats$PropFlow <- flow$PropProtoOmi[match(stats$SubjTime, flow$SubjTime)]

#plot the correlation
allColors <- c("Omicron" = "#7C1D6f", 
               "Omicron And Prototype" = "#DC3977",
               "Prototype" = "#045275")

ggplot(stats, aes(x=Proportion, y = PropFlow))+
  geom_point(size = 2, shape = 21, aes(fill = Booster))+
  geom_abline(intercept=0, slope=1)+
  scale_fill_manual(values = allColors)+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  ylab("Proportion Cross-Reactive By Flow")+
  xlab("Proportion Cross-Reactive By CITESeq")+
  ggtitle("Cross-Reactive Response Comparison - 2-8")+
  theme_classic()+
  theme(plot.title = element_text(size=8), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "06_repertoire_analysis", "tinkering", "CITESeq_Vs_Flow_RemovingIndeterminate.png"),width = 3, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#####
#Plot citeseq signal for postbac poster day
ggplot(df[df$Subject == "4953494948",], aes(x = Proto.RBD.PE, y= XBB.RBD.no.fluor))+
  geom_point(shape = 21, size =1.2, stroke = 0.6, aes(fill = ProtoOmi))+
  scale_x_continuous(trans = "log10", limits = c(1, 1000))+
  scale_y_continuous(trans = "log10", limits = c(1, 1000))+
  ylab("Omicron RBD Binding Signal")+
  xlab("Prototype RBD Binding Signal")+
  labs(fill = "Probe-Binding Label")+
  theme_classic()+
  theme(axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=9,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=9),
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.text = element_text(size = 10),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 10, face="bold"),
        legend.margin=margin(0,0,0,0))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Figure4_ProbeSignalForPostbacPosterDay.png"),width = 4, height = 3, units = "in", device = "png", dpi = 600)
dev.off()
