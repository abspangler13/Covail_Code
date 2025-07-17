#load dependencies
library(Seurat)
library(ggplot2)
library(dplyr)
library(Peptides)
library(here)
library(tidyverse)
library(RColorBrewer)
library(readxl)
library(rstatix)
library(writexl)

#load the flow data
flow <- read_xls(here::here("01_raw-data", "FlowData", "ProtoxBAxBivalentdataUnblindedFM.xls"), sheet= "With Pfizer vaccination")
#infected <- unique(flow$`Subject ID`[flow$infect_flag == "Y"])
#flow <- flow[!(flow$`Subject ID` %in% infected) & flow$Stage == 1,] #as of 240129, this dataset is all already stage 1; this is the flow file with Pfizer vaccination included that Flavio made

#we need to temporarily remove some of the serology-cytometry demultiplexed data bc the SN mappings are incorrect
flow <- flow[!(flow$`Subject ID` %in% c(4952494948, 4952564848)),]

#change the treatment labels to be the same- we can keep a placeholder variable for convenience that will reflect Pfizer vaccination
flow$TrueTreatment <- flow$Treatment
flow <- flow %>%
          mutate(Treatment = case_when(Treatment == "Omicron (Pfizer 1)" ~ "1 Dose Omicron (Moderna)",
                                       Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "1 Dose Omicron + Prototype (Moderna)",
                                       Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "1 Dose Prototype (Moderna)",
                                       TRUE ~ Treatment),
                 Treatment = str_remove(Treatment, " \\(Moderna\\)"),
                 timepoint = case_when(`Time point Guess` == 1 ~ "Day 1",
                                       `Time point Guess` == 15 ~ "Day 15",
                                       `Time point Guess` == 90 ~ "Day 91",
                                       `Time point Guess` == 180 ~ "Day 181"))

# > table(flow$Treatment)
# 1 Dose Omicron 1 Dose Omicron + Prototype           1 Dose Prototype   #looking good!
# 76                         82                         75 
  
# > table(flow$`Time point Guess`, flow$timepoint)
# Day 1 Day 15 Day 181 Day 91
# 1      57      0       0      0
# 15      0     59       0      0
# 90      0      0       0     58
# 180     0      0      59      0. #this also checks out!

#load the citeseq data
df <- read.csv(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis","COVAIL_Metadata_ClusteredAndCleaned_Uninfected_AtLeastSinglePositive_NoNaiveCells.csv"))

#remake plots looking at different RBD+ compartments
#####
#scale variables of interest
#####
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

##R is parsing Pfizer frequencies as non-percentages and not dividing them by 100 by default- let's divide the calculated rows by 100 to correct
#Sum each of the totals:
flow$TotalRBD <- rowSums(flow[,c(56:63, 65, 66, 68, 70)])
flow$TotalProtoOnly <-  rowSums(flow[,c(58, 59, 60, 63)])
flow$TotalBA1Only <- rowSums(flow[,c(65, 66, 68, 70)])
flow$TotalCrossReactive <- rowSums(flow[, c(56, 57, 61, 62)])

#make boxplots
#Create a total RBD boxplot
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalRBD.pdf"))
ggplot(flow, aes(x=timepoint, y=TotalRBD, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Total RBD+ Percentage of IgG+")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,13)+
  ggtitle("Total RBD+ Percentage Over Time By Flow")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=10), 
        axis.title.y = element_text(size=10),
        axis.text.x = element_text(size=12,angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=10))
dev.off()

#Create a boxplot for total BA1+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalBA1Only.pdf"))
ggplot(flow, aes(x=timepoint, y=TotalBA1Only, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Total BA1-Only+ Percentage of IgG+")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,8.1)+
  ggtitle("Total BA1-Only+ Counts Over Time By Flow")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
       plot.title = element_text(size=16), 
       axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=16,angle = 90, hjust=1, vjust=0.5),
       axis.text.y = element_text(size=16))
dev.off()

#Create a boxplot for total Prototype+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalProtoOnly.pdf"))
ggplot(flow, aes(x=timepoint, y=TotalProtoOnly, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Total Proto-Only+ Percentage of IgG+")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,8.1)+
  ggtitle("Total Prototype-Only+ Counts Over Time By Flow")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
       plot.title = element_text(size=12), 
       axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=16, angle = 90, hjust=1, vjust=0.5),
       axis.text.y = element_text(size=16))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalProtoOnly_scalesbetter.pdf"))
ggplot(flow, aes(x=timepoint, y=TotalProtoOnly, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Total Proto-Only+ Percentage of IgG+")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,5)+
  ggtitle("Total Prototype-Only+ Counts Over Time By Flow")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=12), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=16))
dev.off()

#Total Cross reactive response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalCrossReactiveOnly.pdf"))
ggplot(flow, aes(x=timepoint, y=TotalCrossReactive, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Total Cross-Reactive Percentage of IgG+")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,8.1)+
  ggtitle("Total Cross-Reactive Counts Over Time By Flow")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
       plot.title = element_text(size=12), 
       axis.title.y = element_text(size=16),
       axis.text.x = element_text(size=16, angle = 90, hjust=1, vjust=0.5),
       axis.text.y = element_text(size=16))
dev.off()


###Let's do lineplots now
stats <- flow %>%
          group_by(Treatment, timepoint) %>%
          summarize(meanTotalRBD = mean(TotalRBD),
                    meanTotalProto = mean(TotalProtoOnly),
                    meanTotalBA1 = mean(TotalBA1Only),
                    meanTotalCrossReactive = mean(TotalCrossReactive),
                    sdTotalRBD = sd(TotalRBD),
                    sdProto = sd(TotalProtoOnly),
                    sdBA1 = sd(TotalBA1Only),
                    sdCross = sd(TotalCrossReactive))

stats <- flow %>%
  group_by(Treatment, timepoint) %>%
  summarize(medianTotalRBD = median(TotalRBD),
            medianTotalProto = median(TotalProtoOnly),
            medianTotalBA1 = median(TotalBA1Only),
            medianTotalCrossReactive = median(TotalCrossReactive),
            sdTotalRBD = IQR(TotalRBD),
            sdProto = IQR(TotalProtoOnly),
            sdBA1 = IQR(TotalBA1Only),
            sdCross = IQR(TotalCrossReactive))


#Create a total RBD Lineplot
# pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalRBD.pdf"))
# ggplot(stats, aes(x=timepoint, y=meanTotalRBD))+
#   geom_point(aes(shape = Treatment, color = Treatment), size=2, position = position_dodge(0.6))+
#   geom_line(aes(color = Treatment, group =Treatment), position = position_dodge(0.6))+
#   geom_errorbar(aes(ymax = meanTotalRBD + sdTotalRBD, ymin = meanTotalRBD - sdTotalRBD, width = 0.2, color = Treatment), position = position_dodge(0.6))+
#   xlab("Timepoint Post-Boost")+
#   ylab("Mean Total RBD+ Percentage of IgG+")+
#   scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
#   ggtitle("Total RBD+ Percentage Over Time By Flow")+
#   theme_classic() +
#   theme(legend.key.size = unit(0.2, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12,angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10))
# dev.off()
# 
# ggplot(flow, aes(x=timepoint, y=TotalRBD))+
#   geom_point(aes(shape = Treatment, color = Treatment), size=2)+
#   geom_line(aes(color = Treatment, group =`Subject ID`))+
#   xlab("Timepoint Post-Boost")+
#   ylab("Mean Total RBD+ Percentage of IgG+")+
#   scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
#   ggtitle("Total RBD+ Percentage Over Time By Flow")+
#   theme_classic() +
#   theme(legend.key.size = unit(0.2, 'cm'),
#         plot.title = element_text(size=10), 
#         axis.title.y = element_text(size=10),
#         axis.text.x = element_text(size=12,angle = 90, hjust=1, vjust=0.5),
#         axis.text.y = element_text(size=10))
# dev.off()


######
#do statistics now- we'll want to compare proportions of each compartment
#this is a tricky statistical situation- we will want to start with a Kruskal Wallis test to compare the groups at every timepoint
#if we find a significant outcome, we will then want to do a pairwise wilcoxon test as a post-hoc analysis
# kruskalwallis <- data.frame()
# for(i in unique(flow$timepoint)){
#   x <- data.frame(Timepoint = i)
#   
#   x$p.value_TotalRBD <-  kruskal.test(TotalRBD ~ Treatment, data = flow[flow$timepoint == i,])$p.value
#   x$p.value_TotalProtoOnly <-  kruskal.test(TotalProtoOnly ~ Treatment, data = flow[flow$timepoint == i,])$p.value
#   x$p.value_TotalBA1Only <-  kruskal.test(TotalBA1Only ~ Treatment, data = flow[flow$timepoint == i,])$p.value
#   x$p.value_TotalCrossReactive <-  kruskal.test(TotalCrossReactive ~ Treatment, data = flow[flow$timepoint == i,])$p.value
#   
#   kruskalwallis <- bind_rows(x, kruskalwallis)
#   
#   if(nrow(kruskalwallis == length(unique(flow$timepoint)))){
#     kruskalwallis$p.value_totalRBD.adj <- p.adjust(kruskalwallis$p.value_TotalRBD, method="holm")
#     kruskalwallis$p.value_totalProto.adj <- p.adjust(kruskalwallis$p.value_TotalProtoOnly, method="holm")
#     kruskalwallis$p.value_totalBA1.adj <- p.adjust(kruskalwallis$p.value_TotalBA1Only, method="holm")
#     kruskalwallis$p.value_totalCrossReactive.adj <- p.adjust(kruskalwallis$p.value_TotalCrossReactive, method="holm")
#   }
#   
#   rm(x)
# }
# 
# write_xlsx(list(KruskalWallisResults = kruskalwallis, DunnTestResults_Day15 = dunnTest), path = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "KruskalWallisTest_DunnPostHoc_PValues.xlsx"))
# 
# #spoiler alert: none of the adjusted p-values are significant >:( I will do a Dunn's posthoc analysis anyways
# #for Day 15 because the p-values are comparatively small; this is not great practice but might be mildly informative
# dunnTest <- bind_rows(dunn_test(flow[flow$timepoint == "Day 15",], TotalRBD ~ Treatment, p.adjust.method = "bonferroni"),
#                       dunn_test(flow[flow$timepoint == "Day 15",], TotalProtoOnly ~ Treatment, p.adjust.method = "bonferroni"),
#                       dunn_test(flow[flow$timepoint == "Day 15",], TotalCrossReactive ~ Treatment, p.adjust.method = "bonferroni"))
# 
# 
# #save file
# write.csv(placeholder[!is.na(placeholder$Group), colnames(placeholder) != "placeholder"],here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "WilcoxonTests_ComparingFlowPopulations.csv"))
# 
# #it's maybe best to try multiple approaches- to minimize the number of tests we use, let's just use one mixed effects linear model
# test = lm(TotalRBD ~ Treatment, data=flow)

#create a dataset for Flavio to use in Prism
stats <- flow %>%
          select(`Subject ID`, Treatment, `Time point Guess`,TotalRBD, TotalProtoOnly, TotalCrossReactive, TotalBA1Only) %>%
          mutate(`Subject ID` = ifelse(duplicated(paste(`Subject ID`, `Time point Guess`)), paste0(`Subject ID`, "_1"), `Subject ID`)) %>%
          pivot_longer(!c(`Subject ID`, Treatment, `Time point Guess`), names_to = "RBD", values_to= "ProportionOfIGG") %>%
          mutate(RBDTime = paste0(RBD, "_Day_", `Time point Guess`),
                 RBDTime = ifelse(str_detect(RBDTime, "Day_180"), str_replace(RBDTime,"Day_180", "Day__180"), RBDTime)) %>%
          select(`Subject ID`, Treatment ,RBDTime, ProportionOfIGG) %>%
          arrange(RBDTime) %>%
          pivot_wider(names_from =RBDTime, values_from = ProportionOfIGG)

write_xlsx(stats, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "FlowData_ProbeSpecificity_PivotedTable.xlsx"))


#######
#Let's plot fold change
flow <- flow %>%
    arrange(timepoint) %>%
    group_by(`Subject ID`) %>%
    mutate(foldTotalRBD = TotalRBD / TotalRBD[1],
           foldBA1Only = TotalBA1Only / TotalBA1Only[1],
           foldProtoOnly = TotalProtoOnly / TotalProtoOnly[1],
           foldCrossReactive = TotalCrossReactive / TotalCrossReactive[1])

#make boxplots
#Create a total RBD boxplot
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalRBD_FOLDCHANGE.pdf"))
ggplot(flow, aes(x=timepoint, y=foldTotalRBD, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Fold Change Total RBD from Day 1")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ggtitle("Fold Change from Day 1 in Total RBD+")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=10), axis.title.y = element_text(size=10))
dev.off()

#Create a boxplot for total BA1+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalBA1Only_FOLDCHANGE.pdf"))
ggplot(flow, aes(x=timepoint, y=foldBA1Only, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("BA1-Only RBD+ Fold Change")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ggtitle("Fold Change in Total BA1-Only+ Counts")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=10), axis.title.y = element_text(size=10))
dev.off()

#Create a boxplot for total Prototype+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalProtoOnly_FOLDCHANGE.pdf"))
ggplot(flow, aes(x=timepoint, y=foldProtoOnly, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Proto-Only RBD+ Fold Change")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ggtitle("Fold Change in Total Proto-Only+ Counts")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=10), axis.title.y = element_text(size=10))
dev.off()

#Total Cross reactive response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_InfectedRemoved_TotalCrossReactiveOnly_FOLDCHANGE.pdf"))
ggplot(flow, aes(x=timepoint, y=foldCrossReactive, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Cross-Reactive Only RBD+ Fold Change")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ggtitle("Fold Change in Total Cross-Reactive Counts")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=10), axis.title.y = element_text(size=10))
dev.off()

#####
#Let's do a proportion of a proportion of IgG to make an equivalent plot for flow data that we make for CITESeq
stats <- flow %>%
          group_by(Treatment ,`Subject ID`, timepoint) %>%
          summarise(PropBA1 = TotalBA1Only / TotalRBD,
                    PropProto = TotalProtoOnly / TotalRBD,
                    PropCross = TotalCrossReactive / TotalRBD)
          
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_Proportion_InfectedRemoved_CrossReactive.pdf"), width= 10, height = 9)
ggplot(stats, aes(x=timepoint, y=PropCross, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Cross-Reactive Proportion of All RBD+ Cells")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,1)+
  ggtitle("Cross-Reactive RBD+ Proportion Over Time By Flow")+
  facet_grid(cols = vars(Treatment))+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=16), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=18, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=18),
        panel.spacing = unit(2, "lines")
        )+
  guides(fill="none")
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_Proportion_InfectedRemoved_Prototype.pdf"), width= 10, height = 9)
ggplot(stats, aes(x=timepoint, y=PropProto, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Prototype-Only Proportion of All RBD+ Cells")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,1)+
  ggtitle("Prototype-Only RBD+ Proportion Over Time By Flow")+
  facet_grid(cols = vars(Treatment))+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=16), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=18, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=18),
        panel.spacing = unit(2, "lines"))+
  guides(fill="none")
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "FlowData_Proportion_InfectedRemoved_BA1.pdf"), width= 10, height = 9)
ggplot(stats, aes(x=timepoint, y=PropBA1, fill=Treatment))+
  geom_boxplot() +
  xlab("Timepoint Post-Boost")+
  ylab("Prototype-Only Proportion of All RBD+ Cells")+
  scale_x_discrete(limits=c("Day 1", "Day 15", "Day 91", "Day 181"))+
  ylim(0,1)+
  ggtitle("Prototype-Only RBD+ Proportion Over Time By Flow")+
  facet_grid(cols = vars(Treatment))+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=16), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=18, angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=18),
        panel.spacing = unit(2, "lines"))+
  guides(fill="none")
dev.off()

#######
#######
#Now let's work with the CITESeq data- we'll want to make a similar figure

#######
stats <- df %>%
          group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
          summarize(n = n()) %>%
          mutate(Proportion = n / sum(n),
                 Subject = as.character(Subject))

#make boxplots
#Create a boxplot for total BA1+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalBA1Only.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto-Omi+",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width= 0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Omicron+ Among Sequenced Cells Per Donor")+
  ylim(0, 1)+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Omicron-Only+ Sequences Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.8, 'cm'), plot.title = element_text(size=14), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 14, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 14),
        panel.spacing = unit(5, "lines"))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

#Create a boxplot for total Prototype+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalProtoOnly.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Proto-Only+ of Sequences Per Donor")+
  ylim(0, 1)+
  facet_grid(cols=vars(Booster))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Prototype-Only+ Counts Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=15), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 16,),
        panel.spacing = unit(1.5, "lines"))
dev.off()

#Total Cross reactive response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalCrossReactiveOnly.pdf"), width= 10, height = 9)
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster))+
  xlab("Timepoint Post-Boost")+
  ylab("Proportion Cross-Reactive of Sequences Per Donor")+
  ylim(0, 1)+
  facet_grid(cols=vars(Booster))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Proportion Cross-Reactive Cells Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=15), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size = 16, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 16),
        panel.spacing = unit(1.5, "lines"))
dev.off()

#write an excel file for Flavio
s <- stats %>%
          mutate(RBDTime = paste0(adj.ProtoOmi, "_", Timepoint),
                 RBDTime = ifelse(str_detect(RBDTime, "Day 180"), str_replace(RBDTime,"Day 180", "Day_180"), RBDTime)) %>%
          ungroup() %>%
          select(Subject, Booster, RBDTime, Proportion) %>%
          arrange(RBDTime) %>%
          pivot_wider(names_from = RBDTime, values_from = Proportion)
          
write_xlsx(s, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "CITESeq_Probe_Specificity_Proportions.xlsx"))

#####
#Let's do fold change on the Citeseq data
stats <- df %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject))%>%
  ungroup()%>%
  arrange(Timepoint) %>%
  group_by(Subject, adj.ProtoOmi) %>%
  mutate(foldChange = Proportion / Proportion[1])

#make boxplots
#Create a boxplot for total BA1+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalBA1Only_FOLDCHANGE.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto-Omi+",], aes(x=Timepoint, y=foldChange, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width=0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Fold Change in Proportion Omicron+ Among Sequenced Cells Per Donor")+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Fold Change in Proportion Omicron-Only+ Sequences Over Time From Day 1")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=15), axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 16, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 16))
dev.off()

#Create a boxplot for total Prototype+ response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalProtoOnly_FOLDCHANGE.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=foldChange, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width=0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Fold Change in Proportion Proto-Only+ of Sequences Per Donor")+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Fold Change in Proportion Prototype-Only+ Counts Over Time From Day 1")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=15), axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 16, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 16))
dev.off()

#Total Cross reactive response
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CITESeqData_InfectedRemoved_TotalCrossReactiveOnly_FOLDCHANGE.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=foldChange, fill=Booster))+
  geom_boxplot() +
  geom_point(shape=21, aes(fill=Booster), position=position_dodge(width=0.75))+
  xlab("Timepoint Post-Boost")+
  ylab("Fold Change in Proportion Cross-Reactive of Sequences Per Donor From Day 1")+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ggtitle("Fold Change in Proportion Cross-Reactive Counts Over Time")+
  theme_classic() +
  theme(legend.key.size = unit(0.2, 'cm'), plot.title = element_text(size=15), axis.title.y = element_text(size=10),
        axis.text.x = element_text(size = 16, angle = 90, hjust= 1, vjust=0.5),
        axis.text.y = element_text(size = 16))
dev.off()

#####
#Let's look at how clonality and probe specificity correlate- 
nonsinglets <- df$clone_subject_id[duplicated(df$clone_subject_id)]
df$CloneStatus <- ifelse(df$clone_subject_id %in% nonsinglets, "Clonal", "Singlet")

stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>%
          group_by(Booster, Subject, CloneStatus, Timepoint, adj.ProtoOmi) %>%
          summarize(n = n()) %>%
          mutate(Proportion = n / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CrossReactivityAtDay15_ClonesVsSinglets_InfectedRemoved.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+" & stats$Timepoint == "Day 15",], aes(x = Booster, y = Proportion, fill=CloneStatus))+
  geom_boxplot()+
  geom_point(aes(fill = CloneStatus), shape=21, position = position_dodge(0.75))+
  ylab("Proportion of Cells That Are Cross-Reactive")+
  ylim(0,1)+
  ggtitle("Proportion of Cells That Are Cross-Reactive at Day 15")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "PrototypeSpecificAtDay15_ClonesVsSinglets_InfectedRemoved.pdf"))
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-" & stats$Timepoint == "Day 15",], aes(x = Booster, y = Proportion, fill=CloneStatus))+
  geom_boxplot()+
  geom_point(aes(fill = CloneStatus), shape=21, position = position_dodge(0.75))+
  ylab("Proportion of Cells That Are Cross-Reactive")+
  ylim(0,1)+
  ggtitle("Proportion of Cells That Are Cross-Reactive at Day 15")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust =0.5))
dev.off()

#make donut plot thingy - we'll need to do a cumulative summary of each clonal group for each timepoint
df$CloneStatus <- ifelse(df$clone_subject_id %in% nonsinglets, df$clone_subject_id, "Singlet")
crossclones <- df$clone_subject_id[df$adj.ProtoOmi == "Proto+Omi+"]
protoclones <- df$clone_subject_id[df$adj.ProtoOmi == "Proto+Omi-"]

stats <- df[df$adj.ProtoOmi != "Proto-Omi+",] %>% #we need to set the color scheme variable
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

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "ProbeSpecifictyAndClonality_Pooled_BoostersOverTime.pdf"), height=10, width=12)
p <- ggplot(stats)+
  geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=2.9, fill=adj.CloneStatus))+
  coord_polar(theta="y")+
  xlim(c(2,4))+
  scale_fill_manual(values = c("Singlet" = "gray80",
                               "Prototype-Specific" = "skyblue",
                               "Cross-Reactive" = "blue"))+
  ggtitle("Clonality and Probe Specificity of Booster Groups Over Time")+
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

#####
#Let's compare SHM from cross-reactive cells present at timepoints 0 and 15 and see if they mature by day 180
stats <- df %>%
          group_by(Booster, clone_subject_id, Timepoint) %>%
          summarize(n = n()) %>%
          mutate(SimpleTimepoint = ifelse(Timepoint %in% c("Day 0", "Day 15"), "Early", Timepoint),
                 LastingClone = ifelse(length(intersect(c("Early", "Day 180"), unique(SimpleTimepoint))) < 2, "No", "Yes"))
earlylateclones <- unique(stats$clone_subject_id[stats$LastingClone == "Yes"])

df$SimpleTimepoint <- ifelse(df$Timepoint %in% c("Day 0", "Day 15"), "Early", df$Timepoint)

stats <- df[df$clone_subject_id %in% earlylateclones & df$Timepoint != "Day 90" & df$adj.ProtoOmi == "Proto+Omi+",] %>% #let's save it for day 180 only for now- it would be a better timepoint to spot maturation
          group_by(Booster, clone_subject_id, SimpleTimepoint) %>%
          summarise(meanSHM = mean(mu_freq))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "CrossReactiveClones_EarlyvsDay180SHM.pdf"), height=10, width=11)
ggplot(stats[stats$meanSHM < 0.18,], aes(x = SimpleTimepoint, y= meanSHM))+ #removed an outlier from Omicron group
  geom_line(aes(group = clone_subject_id))+
  geom_point(shape =21, aes(fill = Booster))+
  ylab("Mean SHM Per Clonal Group")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("Early", "Day 180"))+
  facet_grid(~Booster)+
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

#also vh genes-
list <- fct_rev(reorder(df$v_call[df$adj.ProtoOmi == "Proto+Omi+"], df$v_call[df$adj.ProtoOmi == "Proto+Omi+"], FUN = length))
omi <- as.character(unique(df$Subject[df$Booster == "Omicron"]))
protoomi <- as.character(unique(df$Subject[df$Booster == "Omicron And Prototype"]))
proto <- as.character(unique(df$Subject[df$Booster == "Prototype"]))

stats <- df[!df$adj.ProtoOmi %in% c("Proto-Omi-", "Proto-Omi+", "Proto+Omi-") & df$Timepoint %in% c("Day 0", "Day 15"),] %>%
  group_by(Subject, Timepoint, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         IGHV = v_call,
         Subject = as.character(Subject)) %>%
  ungroup() %>%
  select(Subject, Timepoint, IGHV, Proportion) %>%
  complete(Subject, Timepoint, IGHV, fill=list(Proportion=0))%>%
  mutate(Booster = case_when(Subject %in% omi ~ "Omicron",
                             Subject %in% protoomi ~ "Omicron And Prototype",
                             Subject %in% proto ~ "Prototype"))



ggplot(stats[stats$Timepoint == "Day 15",], aes(x=IGHV, y=Subject, alpha=Proportion, fill = Booster))+
  geom_tile(color="white", linewidth=1)+
  scale_alpha(limits= c(0,0.15),range = c(0.1, 1))+ #I had to put the upper limit super low- IGHV1-69 dominates every timepoint for every group, and the rest of the reponse is more diverse
  coord_equal()+
  scale_y_discrete(limits= c(proto, protoomi, omi))+
  scale_x_discrete(limits = levels(list))+
  labs(x="IGHV Gene", y="Timepoint", title="IGHV Gene Usage Over Time, Cross-Reactive Cells")+
  #scale_fill_brewer(palette="Blues")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9), axis.text.y = element_text(size=9))
ggsave(filename = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "VHUsageByGroup_OverTime_Tilemap_CrossReactiveOnly.png"),width = 20, height = 10, units = "in", device = "png", dpi = 600)
dev.off()
#

#####
#Making an example plot of CITESeq data for the DMID meeting :)
#we'll also want RBD- cells so let's load in the seurat object
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones.rds"))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbeSpecificity", "PrototypeVsXBBPlotted_ForDMIDMeeting.pdf"), height=6, width=7)
ggplot(seuObj@meta.data, aes(x= Proto.RBD.PE, y = XBB.RBD.no.fluor))+
  geom_point(shape = 21, aes(fill = adj.ProtoOmi))+
  scale_x_continuous(trans = 'log2')+
  scale_y_continuous(trans = 'log2')+
  scale_fill_manual(values = c("Proto-Omi-" = "gray",
                               "Proto+Omi-" = "#00BA38",
                               "Proto-Omi+" = "#F8766D",
                               "Proto+Omi+" = "#619CFF"))+
  xlab("Prototype Counts")+
  ylab("Omicron Counts")+
  ggtitle("Prototype vs XBB.1.5 Sequencing Signal")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=11), 
        axis.text.y = element_text(size=11),
        legend.key.size = unit(0.2, 'cm'))
dev.off()
#####

#####
###Code graveyard
# placeholder <- data.frame(placeholder = "yes")
# 
# for(i in unique(flow$`Time point Guess`)){ #stats between day 0 and day 180 for each booster group
#   x <- data.frame(Group = i)
#   x$Comparison <- "between groups- O vs P"
#   
#   x$wilcox.unpaired.p.value_TotalRBD <- wilcox.test(flow$TotalRBD[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                          flow$TotalRBD[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                          alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalProtoOnly <- wilcox.test(flow$TotalProtoOnly[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                                     flow$TotalProtoOnly[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                     alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalCrossReactive <- wilcox.test(flow$TotalCrossReactive[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                                           flow$TotalCrossReactive[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                           alternative = "two.sided", paired = FALSE)$p.value
#   
#   placeholder <- bind_rows(x, placeholder)
# }
# 
# for(i in unique(flow$`Time point Guess`)){ #stats between day 0 and day 180 for each booster group
#   x <- data.frame(Group = i)
#   x$Comparison <- "between groups- O vs O+P"
#   
#   x$wilcox.unpaired.p.value_TotalRBD <- wilcox.test(flow$TotalRBD[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                                     flow$TotalRBD[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                     alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalProtoOnly <- wilcox.test(flow$TotalProtoOnly[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                                           flow$TotalProtoOnly[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                           alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalCrossReactive <- wilcox.test(flow$TotalCrossReactive[flow$Treatment == "1 Dose Omicron (Moderna)" & flow$`Time point Guess` == i],
#                                                               flow$TotalCrossReactive[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                               alternative = "two.sided", paired = FALSE)$p.value
#   
#   placeholder <- bind_rows(x, placeholder)
# }
# 
# for(i in unique(flow$`Time point Guess`)){ #stats between day 0 and day 180 for each booster group
#   x <- data.frame(Group = i)
#   x$Comparison <- "between groups- P vs O+P"
#   
#   x$wilcox.unpaired.p.value_TotalRBD <- wilcox.test(flow$TotalRBD[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                     flow$TotalRBD[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                     alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalProtoOnly <- wilcox.test(flow$TotalProtoOnly[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                           flow$TotalProtoOnly[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                           alternative = "two.sided", paired = FALSE)$p.value
#   
#   x$wilcox.unpaired.p.value_TotalCrossReactive <- wilcox.test(flow$TotalCrossReactive[flow$Treatment == "1 Dose Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                               flow$TotalCrossReactive[flow$Treatment == "1 Dose Omicron + Prototype (Moderna)" & flow$`Time point Guess` == i],
#                                                               alternative = "two.sided", paired = FALSE)$p.value
#   
#   placeholder <- bind_rows(x, placeholder)
# }