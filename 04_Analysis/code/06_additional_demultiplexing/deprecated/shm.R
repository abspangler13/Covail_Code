#tinkering with shm
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
library(stringdist)
library(tidyseurat)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))

df <- as.data.frame(seuObj@meta.data) %>% filter(adj.ProtoOmi != "Proto-Omi-" & ClusterLabel != "Naive")
df$VHVLMutCount <- df$mu_count + df$LC_mu_count

#complete() is useful for distributive plots, but I can't add booster label in group_by() or it flips out, so i'll add 
o <- unique(df$Subject[df$Booster == "Omicron"])
op <- unique(df$Subject[df$Booster == "Omicron And Prototype"])
p <- unique(df$Subject[df$Booster == "Prototype"])

df$Timepoint <- factor(df$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
df$InfectionTimepoint <- factor(df$InfectionTimepoint, levels = c("Pre-Infection", "Post-Infection"))

#####
#As a QC to SHM labels, let's look at how SHM differs between naive and non-naive cells
seuObj@meta.data$Naive <- ifelse(seuObj@meta.data$seurat_clusters == 8, "Naive", "Non-Naive")
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","SHMOverall_NaiveVsNonNaiveCells.pdf"), height=10, width=11)
ggplot(seuObj@meta.data, aes(x = Naive, y= mu_freq))+
  geom_violin(aes(fill = Naive))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  xlab("Naive Status")+
  ylab("Mutation Frequency")+
  ggtitle("Comparing Naive vs Non-Naive SHM")+
  scale_fill_brewer(palette = "PRGn")+
  theme_classic()+
  theme(legend.text = element_text(size=12),
        plot.title = element_text(size=12), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

######
#first, let's compare SHM for all cells at each timepoint. Are there differences over time?
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","SHMOverall_AllCells_OverTime_FacetedByInfection.pdf"), height=10, width=11)
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
  facet_grid(cols = vars(Infection))+
  ggtitle("SHM Among All Cells - Uninfected And Infected") +
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

#let's look at counts as well
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","SHMOverall_VHVLMuCounts_AllCells_OverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi != "Proto-Omi+",], aes(x = Timepoint, y= VHVLMutCount))+ #removed an outlier from Omicron group
  geom_jitter(shape =21, aes(fill = Timepoint))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  ylab("VH + VL Mutations")+
  xlab("Timepoint")+
  scale_fill_brewer(breaks = c("Day 0", "Day 15", "Day 90", "Day 180"),palette = 1)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols = vars(Infection))+
  ggtitle("SHM Among All Cells - Uninfected And Infected") +
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

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","SHMOverall_VHMuCounts_AllCells_OverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi != "Proto-Omi+",], aes(x = Timepoint, y= mu_count))+ #removed an outlier from Omicron group
  geom_jitter(shape =21, aes(fill = Timepoint))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  ylab("VH Mutations")+
  xlab("Timepoint")+
  scale_fill_brewer(breaks = c("Day 0", "Day 15", "Day 90", "Day 180"),palette = 1)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols = vars(Infection))+
  ggtitle("SHM Among All Cells - Uninfected And Infected") +
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

#let's get a better understanding of the data's spread- let's plot median per individual over time
stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>%
          group_by(Infection, Booster, Subject, Timepoint) %>%
          summarize(medianSHM = median(mu_freq))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","MedianSHMPerIndividual_OverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(stats, aes(x= Timepoint, y= medianSHM))+
  geom_point(shape =21, aes(fill = Booster))+
  geom_line(aes(group = Subject))+
  facet_grid(rows = vars(Booster),cols = vars(Infection))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ylab("Median Proportion SHM Per Subject")+
  xlab("Timepoint")+
  ggtitle("Median SHM Per Subject - Uninfected And Infected") +
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

#okay, it looks like consistently "maturation" occurs by day 15. This is way too soon to experience proper maturation.
#this is probably selective expansion of cross-reactive clones. Let's see if this changes much over time
nonsinglets <- unique(df$clone_subject_id[duplicated(df$clone_subject_id) | duplicated(df$clone_subject_id, fromLast = TRUE)])

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","OverallSHM_ClonesOnly_OverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi != "Proto-Omi+" & df$clone_subject_id %in% nonsinglets,], aes(x = Timepoint, y= mu_freq))+ #removed an outlier from Omicron group
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
  facet_grid(cols = vars(Infection))+
  ggtitle("SHM Among Clones - Uninfected And Infected") +
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

#let's look at just singlets
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","OverallSHM_SingletsOnly_OverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(df[df$adj.ProtoOmi != "Proto-Omi+" & !(df$clone_subject_id %in% nonsinglets),], aes(x = Timepoint, y= mu_freq))+ #removed an outlier from Omicron group
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
  facet_grid(cols = vars(Infection))+
  ggtitle("SHM Among Singlets - Uninfected And Infected") +
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

#We know that clones don't really mature over time at all, so let's see which are being recalled
#Define low vs high SHM clones present before boosting and after. Low = median SHM is lower than median of all cells.
df$SimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 0"), "Preboost", "Postboost")

stats <- df %>%
  group_by(clone_subject_id, SimpleTimepoint) %>%
  summarize(n = n()) %>%
  mutate(LastingClone = ifelse(length(intersect(c("Preboost", "Postboost"), unique(SimpleTimepoint))) < 2, "No", "Yes"))

df$LastingClone <- stats$LastingClone[match(df$clone_subject_id, stats$clone_subject_id)]

#define median shm
medianSHM <- median(df$mu_freq[df$adj.ProtoOmi != "Proto-Omi+" & df$Timepoint == "Day 0"])

#calculate
calculatedMedian <- df %>%
                    group_by(clone_subject_id) %>%
                    summarize(medSHM = median(mu_freq)) %>%
                    mutate(SHMLevel = ifelse(medSHM >= medianSHM, "High", "Low"))

stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>%
          group_by(LastingClone, Subject , clone_subject_id, SimpleTimepoint)%>%
          summarize(n = n()) %>%
          group_by(Subject, SimpleTimepoint) %>%
          mutate(Proportion = n / sum(n))
stats$SHMLevel <- calculatedMedian$SHMLevel[match(stats$clone_subject_id, calculatedMedian$clone_subject_id)]

#plot proportion over time
ggplot(stats[stats$LastingClone == "Yes",], aes(x = SimpleTimepoint, y = Proportion))+
  geom_point(aes(fill = SHMLevel), shape = 21)+
  geom_line(aes(group = clone_subject_id))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y.., color = SimpleTimepoint),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  xlab("Timepoint")+
  ylab("Proportion of Total Response of Individual")+
  scale_x_discrete(limits= c("Preboost", "Postboost"))+
  ggtitle("Change in Proportion of Clones By SHM Level")+
  facet_grid(cols = vars(SHMLevel))+
  theme_classic()+
  theme()
dev.off()


#first, plot SHM of lasting clones vs clones that don't last
ggplot(df[df$clone_subject_id %in% nonsinglets,], aes(x = LastingClone, y = mu_freq))+
  geom_jitter(shape =21, aes(fill = LastingClone))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..),
               position = position_dodge(width = 1),
               linewidth = 1.5)+
  ylab("Proportion SHM")+
  xlab("Lasting vs Non-Lasting Clones")+
  facet_grid(cols = vars(Timepoint))+
  ggtitle("SHM Among All Cells - Uninfected And Infected") +
  theme_classic()+
  theme()
dev.off()

#it's really hard to tell what the dynamic here is. Let's try doing this with just day 0 and day 15
df$LateTimes <- ifelse(df$Timepoint %in% c("Day 90", "Day 180"), "Late", df$Timepoint)
stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>%
          group_by(clone_subject_id, LateTimes) %>%
          summarize(medSHM = median(mu_freq)) %>%
          mutate(PresentAcrossAllTimepoints = ifelse(length(intersect(c("Day 0", "Day 15", "Late"), unique(LateTimes))) == 3, "Yes", "No"))

stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>%
  group_by(Infection, clone_subject_id, Timepoint) %>%
  summarize(medSHM = median(mu_freq)) %>%
  mutate(PresentAcrossAllTimepoints = ifelse(length(intersect(c("Day 0", "Day 15", "Day 90", "Day 180"), unique(Timepoint))) == 4, "Yes", "No"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis","SHM","MedianSHM_LastingClonesOverTime_FacetedByInfection.pdf"), height=10, width=11)
ggplot(stats[stats$PresentAcrossAllTimepoints == "Yes",], aes(x = Timepoint, y=medSHM))+
  geom_point()+
  geom_line(aes(group = clone_subject_id))+
  stat_summary(fun = mean,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y..), color = "red",
                 position = position_dodge(width = 1),
               linewidth = 1.5)+
  ggtitle("Median SHM of Clonal Groups Over Time")+
  facet_grid(cols = vars(Infection))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
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

#overall, it looks like SHM increases mostly at day 15 when considering all cells. This is something
#observable among singlets and clones individually as well.
#When considering only clones that last throughout the response, there may be a slight increase in SHM over time
#However, whether this is proper maturation is unclear, and it stands to reason that there is probably just selective
#pressure for more mutated MBCs
#####


#####
#let's make pyramid plot style thingy
#making a population pyramid-style SHM plot and faceting by group
df <- df %>% mutate(mu_bin = case_when(mu_freq <= 0.025 ~ "<0.025",
                                       mu_freq <= 0.05 ~ "0.025 - 0.05",
                                       mu_freq <= 0.075 ~ "0.05 - 0.075",
                                       mu_freq <= 0.1 ~ "0.075 - 0.1",
                                       mu_freq <= 0.125 ~ "0.1 - 0.125",
                                       mu_freq <= 0.15 ~ "0.125 - 0.15",
                                       mu_freq <= 0.175 ~ "0.15 - 0.175",
                                       mu_freq <= 0.2 ~ "0.175 - 0.2",
                                       mu_freq > 0.2 ~ ">0.20",
                                       TRUE ~ "Hmmmm"))

stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & !df$Timepoint %in% c("Day 0"),] %>%
  group_by(Booster, Infection, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n)))

muRange <- c(-1*max(abs(stats$n)), max(abs(stats$n)))
mu_range_breaks <- pretty(muRange, n = 7)

#absolute numbers
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM","SHMOverall_DistributionDifferences_PyramidPlot.pdf"), height=10, width=11)
ggplot(stats, aes(x = n, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_fill_brewer(palette  = "BuPu")+
  scale_x_continuous(breaks = mu_range_breaks,
                     labels = abs(mu_range_breaks))+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Excluding Day 0")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

#proportion
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.45, 0.45)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Excluding Day 0")+
  xlab("Proportion of All Cells")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()


######parsing out the infected people is sort of hard- let's focus on just the uninfected and see changes over time
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "N",] %>%
  group_by(Booster, Timepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n)))

muRange <- c(-1*max(abs(stats$n)), max(abs(stats$n)))
mu_range_breaks <- pretty(muRange, n = 7)

#absolute numbers
stats$Timepoint <- factor(stats$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_facetedbytime.pdf"), height=10, width=11)
ggplot(stats, aes(x = n, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  scale_x_continuous(breaks = mu_range_breaks,
                     labels = abs(mu_range_breaks))+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Uninfected Only")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

#proportion
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_facetedbytime.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.5, 0.5)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Uninfected Only")+
  xlab("Proportion of All Cells Per Timepoint And Booster")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()


#######now infected only- let's do over time but also by infection timepoint
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "Y",] %>%
  group_by(Booster, Timepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n)))

muRange <- c(-1*max(abs(stats$n)), max(abs(stats$n)))
mu_range_breaks <- pretty(muRange, n = 7)

#absolute numbers
stats$Timepoint <- factor(stats$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_facetedbytime_infectedonly.pdf"), height=10, width=11)
ggplot(stats, aes(x = n, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  scale_x_continuous(breaks = mu_range_breaks,
                     labels = abs(mu_range_breaks))+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

#proportion
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_facetedbytime_infectedonly.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.5, 0.5)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only")+
  xlab("Proportion of Probe-Spec. Cells Per Timepoint And Booster")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

###now by infection timepoint (i.e. pre vs post)
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "Y" & df$Timepoint != "Day 0",] %>%
  group_by(Booster, InfectionTimepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n)))

muRange <- c(-1*max(abs(stats$n)), max(abs(stats$n)))
mu_range_breaks <- pretty(muRange, n = 7)

#absolute numbers
stats$InfectionTimepoint <- factor(stats$InfectionTimepoint, levels = c("Pre-Infection", "Post-Infection"))
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_facetedbytime_infectedonly_prevspostinf.pdf"), height=10, width=11)
ggplot(stats, aes(x = n, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(InfectionTimepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  scale_x_continuous(breaks = mu_range_breaks,
                     labels = abs(mu_range_breaks))+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

#proportion
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_facetedbytime_infectedonly_prevspostinf.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  geom_col()+
  facet_grid(cols = vars(InfectionTimepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.5, 0.5)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only")+
  xlab("Proportion of Probe-Spec. Cells Per Timepoint And Booster")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

#####on an individual basis
###uninfected
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "N",] %>%
  group_by(Subject, Timepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n))) %>%
  ungroup()%>%
  complete(Subject, Timepoint, adj.ProtoOmi, mu_bin, fill = list(Proportion = 0)) %>%
  mutate(Booster = case_when(Subject %in% o ~ "Omicron", #originally, I added booster in grouping, but complete() is not great when doing this
                             Subject %in% op ~ "Omicron And Prototype",
                             Subject %in% p ~ "Prototype",
                             TRUE ~ "whoops"))
  
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_PerSubject.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  stat_summary(fun = median,
               geom = "bar")+
  geom_point(shape = 21)+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.75, 0.75)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Uninfected Only")+
  xlab("Proportion of Cells Per Timepoint And Individual")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

###infected 
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "Y",] %>%
  group_by(Subject, Timepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n))) %>%
  ungroup()%>%
  complete(Subject, Timepoint, adj.ProtoOmi, mu_bin, fill = list(Proportion = 0)) %>%
  mutate(Booster = case_when(Subject %in% o ~ "Omicron", #originally, I added booster in grouping, but complete() is not great when doing this
                             Subject %in% op ~ "Omicron And Prototype",
                             Subject %in% p ~ "Prototype",
                             TRUE ~ "whoops"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_PerSubject_Infected.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  stat_summary(fun = median,
               geom = "bar")+
  geom_point(shape = 21)+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.75, 0.75)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only")+
  xlab("Proportion of Cells Per Timepoint And Individual")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()

##infected but more demultiplexed
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$Infection == "Y",] %>%
  group_by(Subject, InfectionTimepoint, adj.ProtoOmi, mu_bin)%>%
  summarize(n = n()) %>%
  mutate(n = ifelse(adj.ProtoOmi == "Proto+Omi-", n*-1, n),
         Proportion = ifelse(adj.ProtoOmi == "Proto+Omi-", n / sum(n) * -1, n / sum(n))) %>%
  ungroup()%>%
  complete(Subject, InfectionTimepoint, adj.ProtoOmi, mu_bin, fill = list(Proportion = 0)) %>%
  mutate(Booster = case_when(Subject %in% o ~ "Omicron", #originally, I added booster in grouping, but complete() is not great when doing this
                             Subject %in% op ~ "Omicron And Prototype",
                             Subject %in% p ~ "Prototype",
                             TRUE ~ "whoops"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM", "SHMOverall_DistributionDifferences_PyramidPlot_proportion_PerSubject_Infected_ByInfectTime.pdf"), height=10, width=11)
ggplot(stats, aes(x = Proportion, y = mu_bin, fill = adj.ProtoOmi))+
  stat_summary(fun = median,
               geom = "bar")+
  geom_point(shape = 21)+
  facet_grid(cols = vars(InfectionTimepoint), rows = vars(Booster))+
  scale_fill_brewer(palette  = "BuPu")+
  xlim(-0.75, 0.75)+
  scale_y_discrete(limits = c("<0.025", "0.025 - 0.05", "0.05 - 0.075", "0.075 - 0.1", "0.1 - 0.125", "0.125 - 0.15", "0.15 - 0.175", "0.175 - 0.2", ">0.20"))+
  ggtitle("SHM Distribution - Infected Only - Pre vs Post Infection - All Timepoints")+
  xlab("Proportion of Cells Per Timepoint And Individual")+
  geom_vline(xintercept = 0)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 9, angle = 90, hjust= 1, vjust=0.5))
dev.off()
