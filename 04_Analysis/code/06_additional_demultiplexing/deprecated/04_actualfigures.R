#load dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(Peptides)
library(here)
library(tidyverse)
library(RColorBrewer)
library(rstatix)
library(gridExtra)
library(thematic)
library(stringr)
library(writexl)
library(viridisLite)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))

#set to a df so I don't have to reference metadata every time
df <- as.data.frame(seuObj@meta.data)
df <- df %>% filter(Infection == "N", adj.ProtoOmi != "Proto-Omi-", !seurat_clusters %in% c(5,8))
rownames(df) <- df$CELL
write.csv(df ,file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis","COVAIL_Metadata_ClusteredAndCleaned_Uninfected_AtLeastSinglePositive_PublicClonesLabelled_NoNaive.csv"))

#####
#Isotype switching
#let's first look more broadly at isotypes overall and then see if there's any switching among
#cross-reactive cells specifically
df <- df %>% mutate(SimplifiedIsotype = case_when(c_call %in% c("IGHA1", "IGHA2") ~ "IGHA",
                               c_call %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4") ~ "IGHG",
                               c_call %in% c("IGHD", "IGHM") ~ "IGHM/IGHD"))

# > table(df$SimplifiedIsotype, df$c_call)
#             IGHA1 IGHA2 IGHD IGHG1 IGHG2 IGHG3 IGHG4 IGHM
# IGHA        964   152    0     0     0     0     0    0
# IGHG          0     0    0  7029   809   309   985    0
# IGHM/IGHD     0     0    9     0     0     0     0  213

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CrossReactiveOnly_SimplifiedIsotypeUsageOverTime.pdf"), height=10, width=10)
ggplot(df[!is.na(df$c_call) & df$adj.ProtoOmi == "Proto+Omi+",], aes(x= Timepoint, y=1, fill = SimplifiedIsotype, color = SimplifiedIsotype))+
  geom_bar(position = "fill",stat="identity")+
  ylab("Proportion of Sequences")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("Day 0", "Day 15","Day 90", "Day 180"))+
  facet_grid(~Booster)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        strip.text = element_text(size=12),
        panel.spacing = unit(1.1, "lines"))
dev.off()
#####

#####
#making figures comparing clustering between all 3 groups
df <- merge(df, as.data.frame(seuObj@reductions$harmony.wnn.umap@cell.embeddings), by= "row.names")
df$Timepoint <- factor(df$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
rownames(df) <- df$CELL

pdf(file = here::here("04_Analysis","plots","06_repertoire_analysis","DimPlotOverTime_ByBooster.pdf"), width=8, height = 7)
ggplot(df[df$Infection == "N",], aes(x=harmonywnnUMAP_1, y=harmonywnnUMAP_2))+
  geom_point(aes(color= seurat_clusters, fill=seurat_clusters), size=0.3)+
  ggtitle("UMAP Clustering Over Time")+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  theme_classic()+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()
#####

#####
#Comparing SHM from day 0 to day 180 but for all cells
#using timepoint as a factor is not working for me so I'm gonna have to put back to char :(
df$Timepoint <- as.character(df$Timepoint)
df$SummedMu <- df$LC_mu_count + df$mu_count
df$SimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 0"), "PreBoost", "PostBoost")
df$LessSimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", df$Timepoint)

df$Timepoint <- factor(df$Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180"))
#let's just look at cells overall
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHMOverall_AllCells_Day0vsEverythingElsePooled.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi != "Proto-Omi+",], aes(x = Timepoint, y= mu_freq))+ #removed an outlier from Omicron group
  geom_jitter(shape =21, aes(fill = Timepoint))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  ylab("Proportion SHM")+
  xlab("Timepoint")+
  scale_fill_brewer(breaks = c("Day 0", "Day 15", "Day 90", "Day 180"),palette = 1)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  #facet_grid(cols = vars(Booster))+
  ggtitle("SHM Among All Cells - Uninfected") +
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#let's also compare clones regardless of whether or not they end up in both early vs day 180 timepoints
nonsinglets <- df$clone_subject_id[duplicated(df$clone_subject_id)]

stats <- df[df$clone_subject_id %in% nonsinglets,] %>%
  group_by(Booster, clone_subject_id, SimpleTimepoint, adj.ProtoOmi) %>%
  summarize(meanSHM = mean(mu_freq))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM_AllClonesNotJustPersistent_Day0vsLater.pdf"), height=10, width=11)
ggplot(stats[stats$adj.ProtoOmi != "Proto-Omi+",], aes(x = SimpleTimepoint, y= meanSHM))+
  geom_jitter(shape =21, aes(fill = Booster))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  ylab("Mean SHM Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("PreBoost", "PostBoost"))+
  facet_grid(cols = vars(Booster))+
  ggtitle("SHM Among All Clones") +
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

###clonal groups 
stats <- df %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n()) %>%
  mutate(SimpleTimepoint = ifelse(Timepoint %in% c("Day 0"), "Early", "Late"),
         LastingClone = ifelse(length(intersect(c("Early", "Late"), unique(SimpleTimepoint))) < 2, "No", "Yes"))

earlylateclones <- unique(stats$clone_subject_id[stats$LastingClone == "Yes"])

df$SimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 0"), "PreBoost", "PostBoost")

stats <- df[df$clone_subject_id %in% earlylateclones,] %>% #let's save it for day 180 only for now- it would be a better timepoint to spot maturation
  group_by(Booster, adj.ProtoOmi, clone_subject_id, SimpleTimepoint) %>%
  summarise(meanSHM = mean(mu_freq),
            n = n()) %>%
  group_by(clone_subject_id) %>%
  mutate(sum = sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHMOverallByClonalGroups_ClonesPresentDay0AndLater_Day0VsLater.pdf"), height=10, width=11)
ggplot(stats[stats$adj.ProtoOmi != "Proto-Omi+",], aes(x = SimpleTimepoint, y= meanSHM))+ #removed an outlier from Omicron group
  geom_line(aes(group = clone_subject_id))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y.., color = Booster),
               position = position_dodge(width = 1))+
  geom_point(shape =21, aes(fill = Booster))+
  ylab("Mean SHM Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("PreBoost", "PostBoost"))+
  facet_grid(cols = vars(Booster))+
  ggtitle("SHM Among All Persistent Clones") +
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#let's look at summed VH+VL counts
stats <- df[df$clone_subject_id %in% earlylateclones,] %>% #let's save it for day 180 only for now- it would be a better timepoint to spot maturation
  group_by(Booster, adj.ProtoOmi, clone_subject_id, SimpleTimepoint) %>%
  summarise(meanSHM = mean(SummedMu),
            n = n()) %>%
  group_by(clone_subject_id) %>%
  mutate(sum = sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHMOverallByClonalGroups_EarlyVSDay180_SummedCounts.pdf"), height=10, width=11)
ggplot(stats[stats$adj.ProtoOmi != "Proto-Omi+",], aes(x = SimpleTimepoint, y= meanSHM))+ #removed an outlier from Omicron group
  geom_line(aes(group = clone_subject_id))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y.., color = Booster),
               position = position_dodge(width = 1))+
  geom_point(shape =21, aes(fill = Booster))+
  ylab("Mean VH/VL Counts Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("PreBoost", "PostBoost"))+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
  ggtitle("SHM Among Proto+/Cross-Reactive Clones") +
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#again
stats <- df[df$clone_subject_id %in% nonsinglets,] %>%
  group_by(Booster, clone_subject_id, SimpleTimepoint, adj.ProtoOmi) %>%
  summarize(meanSHM = mean(SummedMu))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHMOverall_AllClonesNotJustPersistent_EarlyVSDay180_SummedCounts.pdf"), height=10, width=11)
ggplot(stats[stats$adj.ProtoOmi != "Proto-Omi+",], aes(x = SimpleTimepoint, y= meanSHM))+
  geom_jitter(shape =21, aes(fill = Booster))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1))+
  ylab("Mean VH/VL Counts Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("PreBoost", "PostBoost"))+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
  ggtitle("SHM Among Cross-Reactive Clones") +
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#let's see how this correlates with affinity
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHMvsXBBSignal_Uninfected_CrossReactiveOnly.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi == "Proto+Omi+" & df$Timepoint != "Day 0",], aes(x= Proto.RBD.PE, y=XBB.RBD.no.fluor))+
  geom_point(shape =21, aes(fill = Booster, size = mu_freq, alpha=mu_freq))+
  ggtitle("SHM vs XBB MSI- Uninfected Donors, Cross-Reactive Cells")+
  scale_y_continuous(trans = "log2")+
  scale_x_continuous(trans = "log2")+
  facet_grid(cols = vars(Booster))+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#####

#####
#let's look at public clonotypes
#if we are proposing that the bivalent booster selectively elicits more cross-reactive clones, then the 
#clear implication here is that there is greater selective pressure for cross-reactivity
#since public clonotypes are products of convergent evolution, we'd expect that pressure to select for cross-reactivity
#would increase the relative proportion of cross-reactive clonotypes in the bivalent group

#Here are some notes on these clonotypes that we'll want to look at:
#I guided some of my choices by eye and by this article: https://www.cell.com/immunity/pdf/S1074-7613(22)00142-X.pdf
#this article profiles a lot of clonotypes- it especially convinced me to include 38057_2811 as it's a bit small but matches very well with that the literature reports
crossReactiveCloneIDs <- c("27672_358", "31072_1890", "38057_2811", "7512_1915") #7512 and 27672 are very similar
protoReactiveCloneIDs <- c("32157_635")

#let's make dotplots that show the proportion per group of each of these clonotypes
op <- unique(df$Subject[df$Booster == "Omicron And Prototype"])
p <- unique(df$Subject[df$Booster == "Prototype"])
o <- unique(df$Subject[df$Booster == "Omicron"])

clonotypesOfInterest <- df[df$Timepoint != "Day 0",] %>%
                        group_by(Booster, Subject, PooledCloneID) %>%
                        summarize(n = n()) %>%
                        mutate(Proportion = n / sum(n)) %>%
                        ungroup()%>%
                        select(Subject, PooledCloneID,  Proportion) %>%
                        complete(Subject, PooledCloneID, fill=list(Proportion=0))%>%
                        mutate(Booster = case_when(Subject %in% op ~ "Omicron And Prototype",
                                                   Subject %in% p ~ "Prototype",
                                                   Subject %in% o ~ "Omicron",
                                                   TRUE ~ "zoinks")) %>%
                        filter(PooledCloneID %in% c(crossReactiveCloneIDs, protoReactiveCloneIDs))%>%
                        mutate(Spec = case_when(PooledCloneID %in% crossReactiveCloneIDs ~ "Cross-Reactive",
                                                PooledCloneID %in% protoReactiveCloneIDs ~ "Prototype-Specific"))

#plot
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "PublicClonotypes_ProportionsByVaccinationGroup.pdf"))
ggplot(clonotypesOfInterest, aes(x = PooledCloneID, y= Proportion))+
  geom_boxplot(aes(fill = Booster))+
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  facet_grid(~Spec, scales = "free_x")+
  ggtitle("Public Clonotypes Per Group - Removed Day 0")+
  ylab("Proportion of Total Response Per Subject")+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 8, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(1.1, "lines"),
        strip.text = element_text(size =10))
dev.off()