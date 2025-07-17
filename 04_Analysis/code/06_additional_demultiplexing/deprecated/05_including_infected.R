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

#load in Seurat object
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))

df <- as.data.frame(seuObj@meta.data) %>% filter(adj.ProtoOmi != "Proto-Omi-" & seurat_clusters != 8)

rownames(df) <- df$CELL

#off-topic, but let's make a seurat object that only contains uninfected people for Abby to use
uninfectedSeuObj <- seuObj %>% filter(Infection == "N")
saveRDS(uninfectedSeuObj, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_UninfectedDonors.rds"))

#set subject names so i can use identifiers down the line for faceted graphing purposes
op <- unique(df$Subject[df$Booster == "Omicron And Prototype" & df$Infection == "Y"])
p <- unique(df$Subject[df$Booster == "Prototype" & df$Infection == "Y"])
o <- unique(df$Subject[df$Booster == "Omicron" & df$Infection == "Y"])

#define nonsinglet
nonsinglets <- df$clone_subject_id[duplicated(df$clone_subject_id) | duplicated(df$clone_subject_id, fromLast = T) & df$Infection == "Y"]

df$CloneStatus <- ifelse(df$clone_subject_id %in% nonsinglets, df$clone_subject_id, "Singlet")
crossclones <- df$clone_subject_id[df$adj.ProtoOmi == "Proto+Omi+" & df$Infection == "Y"]
protoclones <- df$clone_subject_id[df$adj.ProtoOmi == "Proto+Omi-" & df$Infection == "Y"]


#let's do some of the basic analyses
#####
#Let's compare probe specificities
stats <- df %>%
  group_by(Booster, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

#first, we'll facet by infection and compare groups
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalBA1Only.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto-Omi+",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Omicron+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Omicron-Only+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalCrossReactiveOnly.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Cross-Reactive Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Cross-Reactive Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalProtoOnly.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Prototype+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(Infection))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Prototype-Only+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

#now that we've demultiplexed the data, let's only look at infected people and split by when they were infected
#first, we'll facet by infection and compare groups
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedOnly_TotalBA1Only_SplitByInfectionDate.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto-Omi+" & stats$Infection == "Y",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Omicron+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Omicron-Only+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalCrossReactiveOnly_SplitByInfectionDate.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+" & stats$Infection == "Y",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_line(aes(group = Subject))+
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Cross-Reactive Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Cross-Reactive Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalProtoOnly_SplitByInfectionDate.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-" & stats$Infection == "Y",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Prototype+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Prototype-Only+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

#now let's look at whether, at a bulk-level, there are proportional differences in cross-reactivity post- infection for both groups
stats <- df[df$Infection == "Y",] %>%
  group_by(Booster, Infection, InfectionRange, Subject, InfectionTimepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))

#first, we'll facet by infection and compare groups
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ProbeSpecificity", "CITESeqData_InfectedIncluded_TotalCrossReactive_BeforeVsAfterInfection.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=InfectionTimepoint, y=Proportion, fill=Booster))+
  geom_boxplot()+
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Relative to Infection")+
  ylab("Proportion Cross-Reactive+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  facet_grid(cols = vars(Booster))+#, rows = vars(InfectionRange))+
  scale_x_discrete(limits=c("Pre-Infection", "Post-Infection"))+
  ggtitle("Proportion Cross-Reactive+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_text(size = 12, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 12),
        panel.spacing = unit(1, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

#####

#####
#let's compare SHM from early and late timepoints- is there any maturation going on?
stats <- df %>%
  filter(Infection == "Y")%>% #filter out uninfected people
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n()) %>%
  mutate(SimpleTimepoint = ifelse(Timepoint %in% c("Day 0"), "Early", Timepoint),
         LastingClone = ifelse(length(intersect(c("Early", "Day 180"), unique(SimpleTimepoint))) < 2, "No", "Yes"))
earlylateclones <- unique(stats$clone_subject_id[stats$LastingClone == "Yes"])

df$SimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 0"), "Day 0", "Late")

stats <- df[df$clone_subject_id %in% earlylateclones & df$Infection == "Y",] %>%
  group_by(Booster, clone_subject_id, SimpleTimepoint) %>%
  summarise(meanSHM = mean(mu_freq))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "AllClones_InfectedIncluded_EarlyvsDay180SHM.pdf"), height=10, width=11)
ggplot(stats, aes(x = SimpleTimepoint, y= meanSHM))+ #removed an outlier from Omicron group
  geom_line(aes(group = clone_subject_id))+
  geom_point(shape =21, aes(fill = Booster))+
  stat_summary(fun = median,
               geom = "errorbar",
               aes(ymax = ..y.., ymin=..y.., color = Booster),
               position = position_dodge(width = 1))+
  ylab("Mean SHM Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("Day 0", "Late"))+
  facet_grid(~Booster)+
  ggtitle("SHM Among Clones - Infected Only") +
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size=10),
        plot.title = element_text(size=20), 
        axis.title.y = element_text(size=16),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size = 15, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 15),
        #panel.spacing = unit(1.5, "lines"),
        strip.text = element_text(size =16))
dev.off()

#####

#####
#let's check what the isotypes look like- we might expect an increase in IgA over time to reflect new mucosal immunity
df <- df %>% mutate(SimplifiedIsotype = case_when(c_call %in% c("IGHA1", "IGHA2") ~ "IGHA",
                                                  c_call %in% c("IGHG1", "IGHG2", "IGHG3", "IGHG4") ~ "IGHG",
                                                  c_call %in% c("IGHD", "IGHM") ~ "IGHM/IGHD"))

# > table(df$SimplifiedIsotype, df$c_call)
#             IGHA1 IGHA2 IGHD IGHG1 IGHG2 IGHG3 IGHG4 IGHM #i've since removed cluster 7 to focus only on memory responses
# IGHA        964   152    0     0     0     0     0    0
# IGHG          0     0    0  7029   809   309   985    0
# IGHM/IGHD     0     0    9     0     0     0     0  213

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "AllCells_SimplifiedIsotypeUsageOverTime_InfectedOnly.pdf"), height=10, width=10)
ggplot(df[!is.na(df$c_call) & df$Infection == "Y",], aes(x= Timepoint, y=1, fill = SimplifiedIsotype, color = SimplifiedIsotype))+
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

#do by individual basis
stats <- df[df$Infection == "Y" & !is.na(df$c_call),] %>%
          group_by(Booster, Subject, Timepoint, SimplifiedIsotype) %>%
          summarize(n =n())%>%
          mutate(Proportion = n / sum(n),
                 IndividualNumber = case_when(Booster == "Omicron" ~ match(Subject, o),
                                              Booster == "Omicron And Prototype" ~ match(Subject, op),
                                              Booster == "Prototype" ~ match(Subject, p)),
                 IndividualNumber = as.character(IndividualNumber))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "IGHAUsageOverTime_PerSubject.pdf"), width=10, height=7)
ggplot(stats[stats$SimplifiedIsotype == "IGHA",], aes(x = Timepoint, y = Proportion, fill = Booster, color = Booster))+
  geom_point(aes(shape = IndividualNumber))+
  geom_line(aes(group = Subject))+
  facet_grid(cols = vars(Booster))+
  theme_linedraw()+
  theme()
dev.off()

#####



#####
#Let's see if there are any striking differences in clonality
#it's time for michel nussenzweig donuts again :)
#I think it'd actually be best to do before and after infection for clonality (rather than day 0 vs 15 vs 90, etc.)
#but we (I) don't know when they were infected, so oh well

#####
#define non-singlets
stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$InfectionRange == "Between Days 15-90",] %>% #we need to set the color scheme variable
  filter(Infection == "Y") %>%
  group_by(Booster, Timepoint, CloneStatus) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n),
         Total = sum(n),
         CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
         adj.CloneStatus = case_when(CloneStatus == "Singlet" ~ "Singlet",
                                     CloneStatus %in% crossclones ~ "Cross-Reactive",
                                     CloneStatus %in% protoclones ~ "Prototype-Specific"),
         adj.CloneStatus = fct(adj.CloneStatus, levels = c("Cross-Reactive", "Prototype-Specific","Singlet")),
         Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180")))%>%
  arrange(adj.CloneStatus)%>%
  mutate(ymax = cumsum(Proportion),
         ymin = c(0, head(ymax, n=-1)))

time <- rep(c("Day 0", "Day 15", "Day 90", "Day 180"), 3)
boost <- c(rep("Omicron",4), rep("Omicron And Prototype", 4), rep("Prototype", 4))
dat_text <- data.frame(Timepoint = time, Booster = boost)

label <- c()
for(j in c("Omicron", "Omicron And Prototype", "Prototype")){
  for(i in c("Day 0", "Day 15", "Day 90", "Day 180")){
    label <- append(label, unique(stats$Total[stats$Booster == j & stats$Timepoint == i]))
  }
}

dat_text$label <- label
dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included","ProbeSpecificity", "ProbeSpecifictyAndClonality_Pooled_BoostersOverTime_InfectedOnlyDays15to90.pdf"), height=10, width=12)
p <- ggplot(stats)+
  geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.9, fill=adj.CloneStatus))+
  coord_polar(theta="y")+
  xlim(c(2,4))+
  scale_fill_manual(values = c("Singlet" = "gray80",
                               "Prototype-Specific" = "skyblue",
                               "Cross-Reactive" = "blue"))+
  ggtitle("Clonality and Probe Specificity Over Time - Infected")+
  guides(fill = "none")+
  facet_grid(cols=vars(Timepoint), rows=vars(Booster))+
  theme_void()+
  theme(plot.title = element_text(hjust=0.5, size = 15),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size =15))+
  geom_text(data = dat_text,
            size = 8,
            mapping = aes(x=-Inf, y=-Inf, label = label),
            hjust = 0.5,
            vjust = 0.5)
print(p)
dev.off()

#let's do this again but also color singlets- bivalent group is doing something interesting
df$CloneStatusSpec <- ifelse(df$CloneStatus == "Singlet", paste0("Singlet_", df$adj.ProtoOmi), df$CloneStatus)

stats <- df[df$adj.ProtoOmi != "Proto-Omi+" & df$InfectionRange == "Between Days 15-90",] %>% #we need to set the color scheme variable
  filter(Infection == "Y") %>%
  group_by(Booster, Timepoint, CloneStatusSpec) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n),
         Total = sum(n),
         CloneStatusSpec = fct_reorder(CloneStatusSpec, Proportion, .desc=TRUE),
         adj.CloneStatus = case_when(CloneStatusSpec == "Singlet_Proto+Omi+" ~ "Cross-Reactive Singlet",
                                     CloneStatusSpec == "Singlet_Proto+Omi-" ~ "Prototype-Specific Singlet",
                                     CloneStatusSpec %in% crossclones ~ "Cross-Reactive Clone",
                                     CloneStatusSpec %in% protoclones ~ "Prototype-Specific Clone"),
         adj.CloneStatus = fct(adj.CloneStatus, levels = c("Cross-Reactive Clone", "Cross-Reactive Singlet","Prototype-Specific Clone", "Prototype-Specific Singlet")),
         Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180")))%>%
  arrange(adj.CloneStatus)%>%
  mutate(ymax = cumsum(Proportion),
         ymin = c(0, head(ymax, n=-1)))

time <- rep(c("Day 0", "Day 15", "Day 90", "Day 180"), 3)
boost <- c(rep("Omicron",4), rep("Omicron And Prototype", 4), rep("Prototype", 4))
dat_text <- data.frame(Timepoint = time, Booster = boost)

label <- c()
for(j in c("Omicron", "Omicron And Prototype", "Prototype")){
  for(i in c("Day 0", "Day 15", "Day 90", "Day 180")){
    label <- append(label, unique(stats$Total[stats$Booster == j & stats$Timepoint == i]))
  }
}

dat_text$label <- label
dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included","ProbeSpecificity", "ProbeSpecifictyAndClonality_Pooled_BoostersOverTime_InfectedOnly_SingletsAlsoLabelledDay15toDay90.pdf"), height=10, width=12)
p <- ggplot(stats)+
  geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.9, fill=adj.CloneStatus))+
  coord_polar(theta="y")+
  xlim(c(2,4))+
  scale_fill_manual(values = c("Cross-Reactive Singlet" = "gray30",
                               "Prototype-Specific Singlet" = "gray70",
                               "Prototype-Specific Clone" = "skyblue",
                               "Cross-Reactive Clone" = "blue"))+
  ggtitle("Clonality and Probe Specificity Over Time - Infected")+
  guides(fill = "none")+
  facet_grid(cols=vars(Timepoint), rows=vars(Booster))+
  theme_void()+
  theme(plot.title = element_text(hjust=0.5, size = 15),
        panel.spacing = unit(0, "lines"),
        strip.text = element_text(size =15))+
  geom_text(data = dat_text,
            size = 8,
            mapping = aes(x=-Inf, y=-Inf, label = label),
            hjust = 0.5,
            vjust = 0.5)
print(p)
dev.off()
#####

#let's make a more tabular version of this so we can see the spread of the data- what proportion of each group is clonal?
stats <- df[df$Infection == "Y",] %>% 
  group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
  summarise(n = n()) %>%
  mutate(MoreThanTwo = ifelse(n >=2, n, 0),
         ProportionClonal = sum(MoreThanTwo) / sum(n),
         IndividualNumber = case_when(Booster == "Omicron" ~ match(Subject, o),
                                      Booster == "Omicron And Prototype" ~ match(Subject, op),
                                      Booster == "Prototype" ~ match(Subject, p)),
         IndividualNumber = as.character(IndividualNumber))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ClonalityProportionByGroup_Overall_InfectedOnly.pdf"), width=10, height=7)
ggplot(stats, aes(x=Timepoint, y=ProportionClonal))+
  stat_summary(fun = mean,
               geom = "bar",
               aes(ymax = ..y.., ymin=..y.., fill = Booster),
               position = position_dodge(width = 1))+
  geom_point(aes(fill = Booster, shape = IndividualNumber))+
  geom_line(aes(group = Subject))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  # scale_fill_manual(values=c(factor1Color, factor2Color))+
  # scale_color_manual(values=c(factor1Color, factor2Color))+
  facet_grid(~ Booster)+
  xlab("Timepoint")+ 
  ylab("Proportion Total RBD+ Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Total RBD+ Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 

#cross-reactive only
stats <- df[df$Infection == "Y" & df$adj.ProtoOmi == "Proto+Omi+",] %>% 
  group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
  summarise(n = n()) %>%
  mutate(MoreThanTwo = ifelse(n >=2, n, 0),
         ProportionClonal = sum(MoreThanTwo) / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ClonalityProportionByGroup_Overall_InfectedOnly_CrossReactiveOnly.pdf"), width=10, height=7)
ggplot(stats, aes(x=Timepoint, y=ProportionClonal))+
  stat_summary(fun = mean,
               geom = "bar",
               aes(ymax = ..y.., ymin=..y.., fill = Booster),
               position = position_dodge(width = 1))+
  geom_point(aes(fill = Booster), shape=21, position= position_dodge(width=1))+
  geom_line(aes(group = Subject))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  # scale_fill_manual(values=c(factor1Color, factor2Color))+
  # scale_color_manual(values=c(factor1Color, factor2Color))+
  facet_grid(~ Booster)+
  xlab("Timepoint")+ 
  ylab("Proportion Cross-Reactive Cells Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Total Cross-Reactive Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 
#

#let's see how well clonality correlates with total cross-reactivity
df$CloneOrNot <- ifelse(df$clone_subject_id %in% nonsinglets, "Clone", "Singlet")

stats <- df[df$Infection == "Y",] %>%
          group_by(Booster, Subject, Timepoint, CloneOrNot) %>%
          summarize(n = n()) %>%
          mutate(ProportionClone = n / sum(n))%>%
          filter(CloneOrNot == "Clone") %>%
          select(!n)

stats2 <- df[df$Infection == "Y",] %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(ProportionCross = n / sum(n))%>%
  filter(adj.ProtoOmi == "Proto+Omi+") %>%
  select(!n)

statsmerged <- bind_cols(stats, stats2)
statsmerged <- statsmerged[,c(1:4,9,5,10)]
colnames(statsmerged) <- c("Booster", "Subject", "Timepoint", "CloneOrNot", "adj.ProtoOmi", "ProportionClone", "ProportionCross")

#plot correlation
statsmerged$Timepoint <-factor(stats$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included", "ClonalityvsCrossReactivityCorrelation_InfectedOnly.pdf"), width=10, height=7)
ggplot(statsmerged, aes(x = ProportionClone, y=ProportionCross, fill = Booster))+
  geom_point(shape = 21)+
  geom_smooth(method='lm', formula= y~x)+
  facet_grid(cols = vars(Timepoint), rows = vars(Booster))+
  theme_linedraw()+
  theme()
dev.off()
#####

#####
#let's look at public clonotypes overall- since it's hard to demultiplex by infection timepoint currently,
#a more broad analysis might be better
crossReactiveCloneIDs <- c("27672_358", "31072_1890", "38057_2811", "7512_1915") #7512 and 27672 are very similar
protoReactiveCloneIDs <- c("32157_635")

#let's make dotplots that show the proportion per group of each of these clonotypes
clonotypesOfInterest <- df[df$Timepoint != "Day 0" & df$Infection == "Y",] %>% #disinclude anything that was pre-existing
  group_by(Booster, Subject, PooledCloneID) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()%>%
  select(Subject, PooledCloneID, Proportion) %>%
  complete(Subject, PooledCloneID, fill=list(Proportion=0))%>%
  mutate(Booster = case_when(Subject %in% op ~ "Omicron And Prototype",
                             Subject %in% p ~ "Prototype",
                             Subject %in% o ~ "Omicron",
                             TRUE ~ "zoinks")) %>%
  filter(PooledCloneID %in% c(crossReactiveCloneIDs, protoReactiveCloneIDs))%>%
  mutate(Spec = case_when(PooledCloneID %in% crossReactiveCloneIDs ~ "Cross-Reactive",
                          PooledCloneID %in% protoReactiveCloneIDs ~ "Prototype-Specific"))

#plot
pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "infected_included","PublicClonotypes_ProportionsByVaccinationGroup_InfectedOnly.pdf"))
ggplot(clonotypesOfInterest, aes(x = PooledCloneID, y= Proportion))+
  geom_boxplot(aes(fill = Booster))+
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  facet_grid(~Spec, scales = "free_x")+
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
#####

#####
#Let's see what the relative composition of the response looks like on a phenotypic basis
# ggplot(df[df$Infection == "Y",], aes(x = Timepoint, y= 1, fill = seurat_clusters, color= seurat_clusters))+
#   geom_bar(stat = "identity", position="fill")+
#   scale_x_discrete(levels = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
#   ylab("Proportion of Sequences")+
#   theme_classic()+
#   theme()
#####

#####
#I'm curious as to how the mean "MFI" of these probes looks- let's plot their distributions to see if there
#are any notable trends- I think it's in our best interest to only look at prototype and XBB; BA.1 signal is *low*
#grab some numerical estimates of the data's spread
statsq <- df %>%
          group_by(Booster, Infection, Timepoint, CloneSubjectIDTimepoint) %>%
          summarize(medMSIProto = quantile(Proto.RBD.PE, 0.70),
                    medMSIXBB = quantile(XBB.RBD.no.fluor, 0.70),
                    medBA1 = quantile(BA1.RBD.PE, 0.70),
                    n=n())

ggplot(df[df$XBB.RBD.no.fluor < 30000,], aes(x = Timepoint, y=XBB.RBD.no.fluor, fill = Booster))+
  geom_boxplot()+
  facet_grid(cols = vars(Infection))+
  theme_classic()+
  theme()
dev.off()
#####

#####
#Grab the clones that meet the following criteria:
#total number of clones over all timepoints: 15
#present at day 180 and one of the previous three timepoints
clones <- df[df$Infection == "N",] %>%
            group_by(Booster, clone_subject_id, Timepoint) %>%
            summarize(n = n()) %>%
            mutate(Total = sum(n),
                   TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & length(intersect(unique(Timepoint), c("Day 0", "Day 15", "Day 90"))) >= 1 ~ "Present",
                                                TRUE ~ "Nope"),
                   Select = Total >= 15 & TimepointCorrect == "Present")

check <- clones[clones$Select == TRUE,] #check outcome
unique(check$clone_subject_id[check$Select == TRUE])
# [1] "5048574848_144_1"    "5048574848_695_111"  "5048574848_906_37"   "4953494948_253_70"  
# [5] "4953494948_323_23"   "4953494948_396_8"    "4953494948_463_56"   "4953494948_53_24"   
# [9] "4953494948_64_29"    "5357484948_1471_110" "5357484948_1499_130" "5357484948_907_299" 
# [13] "5553564848_1640_26"  "4955534848_241_64"   "4955534848_393_108"  "4955534848_483_56"  
# [17] "4955534848_980_84"  