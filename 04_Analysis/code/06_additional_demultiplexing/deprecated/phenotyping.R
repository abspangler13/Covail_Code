#looking at phenotyping over time
#load dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(gridExtra)
library(ggalluvial)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))

df <- as.data.frame(seuObj@meta.data) %>% filter(adj.ProtoOmi != "Proto-Omi-" & ClusterLabel != "Naive") #exclude naive cells
df$Timepoint <- factor(df$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

#######
#kinetics - do we see changes in clusters over time?
stats <- df %>%
  group_by(Infection, Timepoint, ClusterLabel) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

#plot alluvial plot of cluster labels
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_InfectionStatus_ClusterCompositionOverTime.pdf"), width = 10, height =8)
ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Spectral")+
  facet_grid(cols = vars(Infection))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

#add in booster
stats <- df %>%
          group_by(Booster, Infection, Timepoint, ClusterLabel) %>%
          summarize(n= n()) %>%
          mutate(Proportion = n / sum(n))

#plot alluvial plot of cluster labels
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_Uninfected_ClusterCompositionOverTime.pdf"), width = 10, height =8)
ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Spectral")+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

#subdivide by probe specificity
stats <- df[df$Infection != "N",] %>%
  group_by(Booster, adj.ProtoOmi,Timepoint, ClusterLabel) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_Uninfected_ClusterCompositionOverTime_SubdividedByProbeSpecificity.pdf"), width = 10, height =8)
ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Spectral")+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off() #it looks like the IgAs could be weird- a lot of Proto-Omi+ cells remain in this cluster

#compare cross-reactivity among cluster labels
stats <- df %>%
  group_by(Infection, ClusterLabel, adj.ProtoOmi) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))
stats$ClusterLabel <- factor(stats$ClusterLabel, levels = c("Activated IgA Memory","Resting IgA Memory", "AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Resting IgG Memory", "Resting IgG Memory 2", "Unclear"))

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_ProbeSpecificityOfEachCluster.pdf"), width = 10, height =8)
ggplot(stats, aes(y = Proportion, x= ClusterLabel, alluvium = adj.ProtoOmi, fill = adj.ProtoOmi, label = adj.ProtoOmi, stratum = adj.ProtoOmi))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Bu")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(breaks = c("Activated IgA Memory","Resting IgA Memory", "AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Resting IgG Memory", "Resting IgG Memory 2", "Unclear"))+
  facet_grid(cols = vars(Infection))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off() #it looks like the IgAs could be non-specific- a lot of Proto-Omi+ cells remain in this cluster

#signal differences between each cluster
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_ProtoVsXBBByCluster.pdf"), width = 10, height =8)
ggplot(df[df$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = XBB.RBD.no.fluor, fill = ClusterLabel))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_wrap(~ClusterLabel)+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_ProtoVsXBBByCluster_CloneCorrectedCells.pdf"), width = 10, height =8)
ggplot(df[df$ProtoOmi != df$adj.ProtoOmi,], aes(x = Proto.RBD.PE, y = XBB.RBD.no.fluor, fill = adj.ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_wrap(~ClusterLabel)+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

#let's track where they go over time phenotypically- let's select clones that persist over the four timepoints and track their phenotypes
#it's gonna be impossible to do an alluvial plot this way, but what we can do is count a whole clonal group as its most common phenotype at that timepoint
persistent <- df %>%
              group_by(clone_id, Timepoint) %>%
             summarize(n = n()) %>%
              mutate(persistent = ifelse(length(unique(Timepoint)) == 4, TRUE, FALSE)) %>%
              filter(persistent)

stats <- df[df$Infection != "N" & df$clone_id %in% unique(persistent$clone_id),] %>%
  group_by(Timepoint, clone_id,ClusterLabel) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(Timepoint, clone_id, ClusterLabel, fill = list(Proportion = 0))

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_Uninfected_ClonalAlluvial.pdf"), width = 10, height =8)
ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = clone_id, fill = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Spectral")+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()
#alluvial plot above will probably not work as is. The issue is that there are multiple starting phenotypes a clonal group could
#have, and they lead to multiple phenotypes. We need each individual ("alluvium") to have a single defined value (i.e. clone_id).
#One answer could be to make a set of data where the most common phpenotype of a clonalgroup is labelled and then plotted

######
#Probe specificity by timepoint and group
stats <- df %>% filter(Infection == "N") %>%
  group_by(Booster, Timepoint, ClusterLabel, adj.ProtoOmi) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))
stats$ClusterLabel <- factor(stats$ClusterLabel, levels = c("Activated IgA Memory","Resting IgA Memory", "AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Resting IgG Memory", "Resting IgG Memory 2", "Unclear"))

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "Phenotyping_ProbeSpecOverTimePerClusterByBooster.pdf"), width = 10, height =10)
ggplot(stats, aes(y = Proportion, x= ClusterLabel, alluvium = adj.ProtoOmi, fill = adj.ProtoOmi, label = adj.ProtoOmi, stratum = adj.ProtoOmi))+
  geom_flow()+
  geom_stratum()+
  scale_fill_brewer(palette = "Bu")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(breaks = c("Activated IgA Memory","Resting IgA Memory", "AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Resting IgG Memory", "Resting IgG Memory 2", "Unclear"))+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

###why do we see *any* activation of prototype-specific cells in the omicron booster?
checkProbes <- df %>% filter(Booster == "Omicron" & ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"))
checkProbes$Mismatch <- checkProbes$adj.ProtoOmi != checkProbes$ProtoOmi

#plot the signals of these 
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "OmicronBooster_AM1AM2_ProtoVsXBB.pdf"), width = 14, height =4)
ggplot(checkProbes[checkProbes$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = XBB.RBD.no.fluor, fill = ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

#plot with BA.1
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "OmicronBooster_AM1AM2_ProtoVsBA1.pdf"), width = 14, height =4)
ggplot(checkProbes[checkProbes$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = BA1.RBD.PE, fill = ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

##what does this look like for prototype vaccination?
#plot with BA.1
checkProbes <- df %>% filter(Booster == "Prototype" & ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"))
checkProbes$Mismatch <- checkProbes$adj.ProtoOmi != checkProbes$ProtoOmi

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "ProtoBooster_AM1AM2_ProtoVsBA1.pdf"), width = 14, height =4)
ggplot(checkProbes[checkProbes$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = BA1.RBD.PE, fill = ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

##plot with adjusted values
#plot with BA.1
checkProbes <- df %>% filter(Booster == "Omicron" & ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"))
checkProbes$Mismatch <- checkProbes$adj.ProtoOmi != checkProbes$ProtoOmi

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "OmicronBooster_AM1AM2_adjProtoVsBA1.pdf"), width = 14, height =4)
ggplot(checkProbes[checkProbes$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = BA1.RBD.PE, fill = ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()

##what does this look like for prototype vaccination?
#plot with BA.1
checkProbes <- df %>% filter(Booster == "Prototype" & ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"))
checkProbes$Mismatch <- checkProbes$adj.ProtoOmi != checkProbes$ProtoOmi

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "phenotyping", "ProtoBooster_AM1AM2_adjProtoVsBA1.pdf"), width = 14, height =4)
ggplot(checkProbes[checkProbes$XBB.RBD.no.fluor < 30000,], aes(x = Proto.RBD.PE, y = BA1.RBD.PE, fill = ProtoOmi))+
  geom_point(shape = 21)+
  scale_fill_brewer(palette = "Spectral")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  theme(axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =12))
dev.off()
