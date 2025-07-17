library(Seurat)
library(sessioninfo)
library(tidyseurat)
library(tidyverse)
library(here)
library(ggplot2)
library(dplyr)
library(writexl)


#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))

#write a nice .csv file with all of the metadata- cleanest version so far
write.csv(seuObj@meta.data, file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis","COVAIL_Metadata_ClusteredAndCleaned.csv"))

#write a seurat object that is totally clean- no naive cells and no more grime (cluster 5)
cleaned <- seuObj %>% filter(ClusterLabel != "Naive" & adj.ProtoOmi != "Proto-Omi-")
saveRDS(cleaned, file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti_NoNaiveCells.rds"))

#set to a df so I don't have to reference metadata every time
metadata <- cleaned@meta.data[cleaned@meta.data$Infection == "N",]

#remove proto-omi-
write.csv(metadata, file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis","COVAIL_Metadata_ClusteredAndCleaned_Uninfected_AtLeastSinglePositive_NoNaiveCells.csv"))

#####
#plot proportions of VH usage over time for each group
#define VH genes that are significant in number:
vhGenes <- c("IGHV1-69", "IGHV3-30", "IGHV3-23", "IGHV4-34", "IGHV5-51","IGHV4-39","IGHV1-46")
metadata$VHGene <- ifelse(metadata$v_call %in% vhGenes, metadata$v_call, "Other")
metadata$VHGene <- reorder(metadata$VHGene, metadata$VHGene, FUN =length)

#plot them as a barplot sort of method
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi != "Proto-Omi-",], aes(fill=VHGene, y=1, x=Booster, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")), switch = "x")+
  ggtitle("Probe-Positive Populations Over Time")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

##
#Let's make a heatmap sort of thing just so that we have it
list <- fct_rev(reorder(metadata$v_call, metadata$v_call, FUN = length))

stats <- metadata[metadata$adj.ProtoOmi != "Proto-Omi-",] %>%
          group_by(Booster, Timepoint, v_call) %>%
          summarize(n = n()) %>%
          mutate(Proportion = n / sum(n),
                 IGHV = v_call,
                 GroupTimepoint = paste0(Booster, " ", Timepoint))

ggplot(stats, aes(x=IGHV, y=GroupTimepoint, alpha=Proportion, fill = Booster))+
  geom_tile(color="white", linewidth=1)+
  scale_alpha(limits= c(0,0.05),range = c(0.2, 1))+ #I had to put the upper limit super low- IGHV1-69 dominates every timepoint for every group, and the rest of the reponse is more diverse
  coord_equal()+
  scale_y_discrete(limits=rev(c("Omicron Day 0", "Omicron Day 15", "Omicron Day 90", "Omicron Day 180",
                              "Omicron And Prototype Day 0", "Omicron And Prototype Day 15", "Omicron And Prototype Day 90", "Omicron And Prototype Day 180",
                              "Prototype Day 0", "Prototype Day 15", "Prototype Day 90", "Prototype Day 180")))+
  scale_x_discrete(limits = levels(list))+
  labs(x="IGHV Gene", y="Timepoint", title="IGHV Gene Usage Over Time")+
  #scale_fill_brewer(palette="Blues")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=9), axis.text.y = element_text(size=9))
ggsave(filename = here::here("04_Analysis", "plots", "06_repertoire_analysis", "VHUsageByGroup_OverTime_Tilemap_comparison.png"),width = 20, height = 10, units = "in", device = "png", dpi = 600)
dev.off()

#####Split VH usage plots by specificity
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_ProtoPositive.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto+Omi-",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ Booster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Proto+ Cells")+
  theme_bw()+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_OmicronPositive.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto-Omi+",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ Booster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Omicron+ Cells")+
  theme_bw()+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_OmicronPrototypePositive.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto+Omi+",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ Booster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Cross-Reactive Cells")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

#do the same but with O+P and O merged
metadata$MergedBooster <- ifelse(metadata$Booster == "Prototype", "Prototype", "Omicron")

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_ProtoPositive_Merged.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto+Omi-",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ MergedBooster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Proto+ Cells")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_OmicronPositive_Merged.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto-Omi+",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ MergedBooster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Omicron+ Cells")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "GroupedStackedBarplot_VHUsageByGroup_OverTime_OmicronPrototypePositive_Merged.pdf"))
ggplot(metadata[metadata$adj.ProtoOmi == "Proto+Omi+",], aes(fill=VHGene, y=1, x=Timepoint, color=VHGene))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of RBD+ Cells")+
  xlab("Timepoint and Booster Group")+
  facet_grid(~ MergedBooster, switch = "x")+ #factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
  ggtitle("Probe-Positive Populations Over Time, Cross-Reactive Cells")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90,  hjust = 0.95, vjust = 0.2), axis.text.y=element_text(size=12),
        strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01, "cm"))
dev.off()
#####

#####
#Look at clonality per individual over time for each group
stats <- metadata[metadata$adj.ProtoOmi != "Proto-Omi-",] %>% 
          group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
          summarise(n = n()) %>%
          mutate(MoreThanTwo = ifelse(n >=2, n, 0),
                 ProportionClonal = sum(MoreThanTwo) / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalityProportionByGroup_Overall.pdf"), width=10, height=7)
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
  ylab("Proportion Total RBD+ Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Total RBD+ Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 

#proto only
stats <- metadata[metadata$adj.ProtoOmi == "Proto+Omi-",] %>% 
  group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
  summarise(n = n()) %>%
  mutate(MoreThanTwo = ifelse(n >=2, n, 0),
         ProportionClonal = sum(MoreThanTwo) / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalityProportionByGroup_ProtoOnly.pdf"), width=10, height=7)
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
  ylab("Proportion Prototype+ Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Prototype+/Omicron- Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 

#omi only
stats <- metadata[metadata$adj.ProtoOmi == "Proto-Omi+",] %>% 
  group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
  summarise(n = n()) %>%
  mutate(MoreThanTwo = ifelse(n >=2, n, 0),
         ProportionClonal = sum(MoreThanTwo) / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalityProportionByGroup_OmicronOnly.pdf"), width=10, height=7)
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
  ylab("Proportion Omicron+ Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Prototype-/Omicron+ Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 

#finally, cross reactive only
stats <- metadata[metadata$adj.ProtoOmi == "Proto+Omi+",] %>% 
  group_by(Booster, Subject, Timepoint, clone_subject_id) %>%
  summarise(n = n()) %>%
  mutate(MoreThanTwo = ifelse(n >=2, n, 0),
         ProportionClonal = sum(MoreThanTwo) / sum(n))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalityProportionByGroup_CrossReactiveOnly.pdf"), width=10, height=7)
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
  ylab("Proportion Cross-Reactive Belonging to a Clonal Group")+
  ggtitle("Clonal Proportion of Cross-Reactive Cells")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7), axis.text.y = element_text(size=9))
dev.off() 

#####
#Let's take a look at the mean SHM divided by vaccination group at day 0 and day 180
#let's first define what is "Activated" and what isn't- what's activated follows two criteria:
#1. at least 5 cells present at day 15
#2. at least a 2-fold increase at day 15
#once you apply these criteria, clones get really low, so I'll define 90 and 180 as a pooled "late" timepoint
metadata$Timepoint.adj <- case_when(metadata$Timepoint %in% c("Day 90", "Day 180") ~ "Day 90/180",
                                    metadata$Timepoint == "Day 15" ~ "Day 15",
                                    metadata$Timepoint == "Day 0" ~ "Day 0")

stats <- metadata %>%
          group_by(clone_subject_id, Booster, Timepoint.adj) %>%
          summarize(meanSHM = mean(mu_freq),
                    n = n()) %>%
          ungroup() %>%
          arrange(Timepoint.adj) %>%
          group_by(clone_subject_id) %>%
          mutate(FoldMeanSHM = meanSHM / meanSHM[1],
                 FoldMeanCounts = n / n[1])
          
large <- stats %>% filter(Timepoint.adj %in% c("Day 0", "Day 90/180")) %>% filter(duplicated(clone_subject_id) | duplicated(clone_subject_id, fromLast =  T))  

#make an overall file for Prism stats
s <- large %>%
      select(clone_subject_id, Booster, Timepoint.adj, meanSHM)%>%
      pivot_wider(names_from = Timepoint.adj, values_from = meanSHM)

write_xlsx(s, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "MeanSHM_StartAndLaterTimepoints_Overall_Pivoted.csv"))

#create a separate dataset with just activated response
activated <- stats %>% 
              filter(clone_subject_id %in% stats$clone_subject_id[(stats$Timepoint.adj == "Day 15") & (stats$n >= 4) & (stats$FoldMeanCounts >= 2)]) %>%
              filter(Timepoint.adj %in% c("Day 0", "Day 90/180")) %>% filter(duplicated(clone_subject_id) | duplicated(clone_subject_id, fromLast =  T))  

#write for an activated-specific file
s <- activated %>%
  select(clone_subject_id, Booster, Timepoint.adj, meanSHM)%>%
  pivot_wider(names_from = Timepoint.adj, values_from = meanSHM)

write_xlsx(s, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "MeanSHM_StartAndLaterTimepoints_ActivatedOnly.xlsx"))

#make plots
pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "MeanSHMPerClonalGroupPresentAtStartAndLaterTimepoints.pdf"), width=7, height=10)
ggplot(large, aes(x=Timepoint.adj, y=meanSHM, group=clone_subject_id, color=Booster))+
  geom_point()+
  geom_line()+
  xlab("Timepoint")+
  ylab("Mean Mu Per Clonal Group")+
  ggtitle("Mean SHM Over Time Per Clonal Group")+
  facet_grid(~ Booster)+
  theme_classic()
dev.off() 

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "MeanSHMPerClonalGroupPresentAtStartAndLaterTimepoints_ActivatedOnly.pdf"), width=7, height=10)
ggplot(activated, aes(x=Timepoint.adj, y=meanSHM, group=clone_subject_id, color=Booster))+
  geom_point()+
  geom_line()+
  xlab("Timepoint")+
  ylab("Mean Mu Per Clonal Group")+
  ggtitle("Mean SHM Over Time Per Clonal Group")+
  facet_grid(~ Booster)+
  theme_classic()
dev.off() 

activated <- activated %>%
          ungroup() %>%
          arrange(Timepoint.adj) %>%
          group_by(clone_subject_id) %>%
          mutate(FoldMeanSHM = meanSHM / meanSHM[1])

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "FoldChangeSHM_PerClonalGroupPresentAtStartAndLaterTimepoints_ActivatedOnly.pdf"), width=7, height=10)
ggplot(activated, aes(x=Timepoint.adj, y=FoldMeanSHM, group=clone_subject_id, color=Booster))+
  geom_point()+
  geom_line()+
  xlab("Timepoint")+
  ylab("Fold Change in Mean Mu Per Clonal Group From Day 0 to Day 180")+
  ggtitle("Fold Change in SHM Over Time")+
  facet_grid(~ Booster)+
  theme_classic()
dev.off()

#we're not seeing much, but this might be because what we see get expanded at day 15 isn't necessarily what we see at day 0. Check between day 15 and day 90/180?
stats <- metadata %>%
  group_by(clone_subject_id, adj.ProtoOmi, Booster, Timepoint.adj) %>%
  summarize(meanSHM = mean(mu_freq),
            n = n()) %>%
  ungroup() %>%
  arrange(Timepoint.adj) %>%
  group_by(clone_subject_id) %>%
  filter(adj.ProtoOmi != "Proto-Omi+")

day15SHM <- stats %>% filter(!(clone_subject_id %in% stats$clone_subject_id[stats$Timepoint.adj == "Day 0"])) %>%
              filter(Timepoint.adj %in% c("Day 15", "Day 90/180")) %>% filter(duplicated(clone_subject_id) | duplicated(clone_subject_id, fromLast =  T))  

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "SHM_PerClonalGroupPresentOnlyAtDay15_SeparatedBySpecificity.pdf"), width=8, height=10)
ggplot(day15SHM, aes(x=Timepoint.adj, y=meanSHM, group=clone_subject_id, color=Booster))+
  geom_point()+
  geom_line()+
  xlab("Timepoint")+
  ylab("Mutation Frequency")+
  ggtitle("Change in SHM for Clones Not Present at Day 0")+
  facet_grid(cols=vars(Booster), rows=vars(adj.ProtoOmi))+
  theme_classic()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))
dev.off()


######
#We'll want to make the michel nussenzweig donuts- we might want to do this iteratively across every donor
#not sure if we really care about pooled clonality
# metadata$SubjGroup <- paste0(metadata$Booster, ", ", metadata$Subject)
# 
# pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "NussenzweigClonalityPies.pdf"))
# for(i in unique(metadata$SubjGroup)){
# placeholder <- metadata[metadata$SubjGroup == i,]
# 
# nonSinglets <- unique(placeholder$clone_subject_id[duplicated(placeholder$clone_subject_id) | duplicated(placeholder$clone_subject_id, fromLast=T)])
# 
# placeholder$CloneStatus <- ifelse(placeholder$clone_subject_id %in% nonSinglets, placeholder$clone_subject_id, "Singlet")          
# 
# placeholder <- placeholder %>%
#   group_by(Timepoint, CloneStatus) %>%
#   summarize(n= n()) %>%
#   mutate(Proportion = n / sum(n),
#          Total = sum(n),
#          CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
#          Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180")),
#          ymax = cumsum(Proportion),
#          ymin = c(0, head(ymax, n=-1)))
# 
# label <- unique(placeholder$Total[order(placeholder$Timepoint)])
# time <- c("Day 0", "Day 15", "Day 90", "Day 180")
# dat_text <- data.frame(label = label, Timepoint = time)
# dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
# 
# p <- ggplot(placeholder)+
#   geom_rect(color= "black", linewidth=0.2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=CloneStatus))+
#   coord_polar(theta="y")+
#   xlim(c(2,4))+
#   scale_fill_manual(values = c("#FFFFFF",viridis(length(unique(placeholder$CloneStatus))-1, begin = 0, end = 1, option = "turbo")))+
#   ggtitle(paste0("Group, Subject: ",i))+
#   guides(fill = "none")+
#   facet_grid(cols=vars(Timepoint))+
#   theme_void()+
#   theme(plot.title = element_text(hjust=0.5))+
#   geom_text(data = dat_text,
#             mapping = aes(x=-Inf, y=-Inf, label = label),
#             hjust = 0.5,
#             vjust = 0.5)
# 
# print(p)
#   }
# dev.off()
# rm(p)
# rm(placeholder)

#now schism the above by reactivity
splitDF <- split(metadata, metadata$adj.ProtoOmi)
splitDF <- splitDF[2:3]
for(j in 1:length(splitDF)){
  
  filename <- paste0(names(splitDF)[j],"_NussenzweigStyleDonuts.pdf")
  df <- splitDF[[j]]
  
  pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", filename))
  for(i in unique(df$Subject)){
    placeholder <- df[df$Subject == i,]
    
    nonSinglets <- unique(placeholder$clone_subject_id[duplicated(placeholder$clone_subject_id) | duplicated(placeholder$clone_subject_id, fromLast=T)])
    placeholder$CloneStatus <- ifelse(placeholder$clone_subject_id %in% nonSinglets, placeholder$clone_subject_id, "Singlet")          
    
    placeholder$Timepoint <- ifelse(placeholder$Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", placeholder$Timepoint)
    
    calcs <- placeholder %>%
              group_by(CloneStatus, Timepoint) %>%
              summarize(n = n()) %>%
              mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                                  length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Day 0 Expanded",
                                  length(unique(Timepoint)) > 1 ~ "Expanded",
                                  TRUE ~ "Single Timepoint"))
    
    placeholder <- placeholder %>%
      group_by(Timepoint, CloneStatus) %>%
      summarize(n= n()) %>%
      mutate(Proportion = n / sum(n),
             Total = sum(n),
             CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
             adj.CloneStatus = case_when( CloneStatus == "Singlet" ~ "Singlet",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Day 0 Expanded"] ~ "Day 0 Expanded",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Expanded"] ~ "Expanded",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Single Timepoint"] ~ "Single Timepoint"),
             adj.CloneStatus = fct(adj.CloneStatus, levels = c("Day 0 Expanded", "Expanded", "Single Timepoint", "Singlet")),
             Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90/180")))%>%
      arrange(adj.CloneStatus)%>%
      mutate(ymax = cumsum(Proportion),
             ymin = c(0, head(ymax, n=-1)))
  
    #placeholder$CloneStatus <- factor(placeholder$CloneStatus, levels= c("Singlet","Expanded", "Day 0 Expanded", "Single Timepoint"))
    label <- c()
    label[1] <- unique(placeholder$Total[placeholder$Timepoint == "Day 0"])
    label[2] <- unique(placeholder$Total[placeholder$Timepoint == "Day 15"])
    label[3] <- unique(placeholder$Total[placeholder$Timepoint == "Day 90/180"])
    time <- c("Day 0", "Day 15", "Day 90/180")
    dat_text <- data.frame(label = label, Timepoint = time)
    dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90/180"))
    

    p <- ggplot(placeholder)+
      geom_rect(color= "black", linewidth=0.2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=adj.CloneStatus))+
      coord_polar(theta="y")+
      xlim(c(2,4))+
      scale_fill_manual(values = c("Singlet" = "#FFFFFF",
                                   "Day 0 Expanded" = "green4",
                                   "Expanded" = "#88CCEE",
                                   "Single Timepoint" = "gray80"))+
      ggtitle(paste0("Group, Subject: ",i, " ", names(splitDF)[j]))+
      guides(fill = "none")+
      facet_grid(cols=vars(Timepoint))+
      theme_void()+
      theme(plot.title = element_text(hjust=0.5))+
      geom_text(data = dat_text,
                mapping = aes(x=-Inf, y=-Inf, label = label),
                hjust = 0.5,
                vjust = 0.5)
    
    print(p)
  }
  dev.off()
  rm(p)
}
rm(splitDF)
rm(p)
rm(calcs)
#####

#####
#Now we need to make a graph in a tabular format- compare timepoints (day 0 vs 15 vs 90/180)
#label clones present at day 0 as one color vs everything new at other times, separate by reactivity
nonSinglets <- unique(metadata$clone_subject_id[duplicated(metadata$clone_subject_id)])
metadata$CloneStatus <- ifelse(metadata$clone_subject_id %in% nonSinglets, metadata$clone_subject_id, "Singlet")          

clones <- metadata[metadata$CloneStatus != "Singlet" & metadata$adj.ProtoOmi != "Proto-Omi+",]
clones$StartingClones <- ifelse(clones$clone_subject_id %in% clones$clone_subject_id[clones$Timepoint == "Day 0"], "Day 0", "Later")
clones$PooledTimepoints <- case_when( clones$Timepoint == "Day 90" ~ "Day 90/180",
                                      clones$Timepoint == "Day 180" ~ "Day 90/180",
                                      TRUE ~ clones$Timepoint)


placeholder <- clones %>%
                group_by(Subject, Booster, adj.ProtoOmi, PooledTimepoints, StartingClones) %>%
                summarize(n = n()) %>%
                mutate(Proportion = n / sum(n),
                       adj.Proportion = ifelse(Booster == "Omicron", Proportion / 4, Proportion / 5 ))

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProportionOfClonesNotFoundAtDay0.pdf"))
ggplot(placeholder[placeholder$StartingClones != "Day 0" & placeholder$PooledTimepoints != "Day 0",], aes(x = PooledTimepoints, y = Proportion, fill = Booster))+
  geom_boxplot()+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
  ylab("Proportion of Clones Not Found at Day 0")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

pdf(here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProportionOfClonesThatAreFoundAtDay0.pdf"))
ggplot(placeholder[placeholder$StartingClones == "Day 0",], aes(x = PooledTimepoints, y = Proportion, fill = Booster))+
  geom_boxplot()+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi))+
  ylab("Proportion of Clones Found at Day 0")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
dev.off()

#Plot specificity changes over time in Clones that overlap with day 0
startClones <- clones[clones$StartingClones == "Day 0",] %>%
                group_by(Booster, Subject, PooledTimepoints, adj.ProtoOmi) %>%
                summarize(n = n()) %>%
                mutate(Proportion = n / sum(n),
                       adj.Proportion = ifelse(Booster == "Omicron", Proportion / 4, Proportion / 5 ))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbePositivePopulations_AmongDay0ClonesOverTime.pdf"),width=8, height = 10)
ggplot(startClones, aes(x = PooledTimepoints, y = adj.Proportion)) +
  geom_bar(stat="identity", aes(y=adj.Proportion, fill=adj.ProtoOmi,  color=adj.ProtoOmi))+
  scale_fill_manual(values = c("Proto+Omi-" = "#D7C49EFF", "Proto+Omi+" = "#343148FF"))+
  scale_color_manual(values = c("Proto+Omi-" = "#D7C49EFF", "Proto+Omi+" = "#343148FF"))+
  facet_grid(cols = vars(Booster))+
  xlab("Timepoint")+
  ylab("Proportion of Probe+ Cells Of Clones Detected Since Day 1")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=16), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16,angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=16))
dev.off()

#clones that don't overlap with day 0
startClones <- clones[clones$StartingClones != "Day 0",] %>%
  group_by(Booster, Subject, PooledTimepoints, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         adj.Proportion = ifelse(Booster == "Omicron", Proportion / 4, Proportion / 5 ))

pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ProbePositivePopulations_AmongDClonesNotPresentAtDay0.pdf"),width=8, height = 10)
ggplot(startClones, aes(x = PooledTimepoints, y = adj.Proportion)) +
  geom_bar(stat="identity", aes(y=adj.Proportion, fill=adj.ProtoOmi,  color=adj.ProtoOmi))+
  scale_fill_manual(values = c("Proto+Omi-" = "#D7C49EFF", "Proto+Omi+" = "#343148FF"))+
  scale_color_manual(values = c("Proto+Omi-" = "#D7C49EFF", "Proto+Omi+" = "#343148FF"))+
  facet_grid(cols = vars(Booster))+
  xlab("Timepoint")+
  ylab("Proportion of Probe+ Cells Of Clones Detected Since Day 1")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=16), 
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16,angle = 90, hjust=1, vjust=0.5),
        axis.text.y = element_text(size=16))
dev.off()

#very quickly make a spreadsheet for Flavio to do stats with
startClones <- clones[clones$StartingClones == "Day 0",] %>%
  group_by(Booster, Subject, PooledTimepoints, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(adj.ProtoOmi == "Proto+Omi+") %>%
  select(Booster, Subject, Subject, PooledTimepoints, Proportion) %>%
  pivot_wider(names_from= PooledTimepoints, values_from=Proportion)

#write xlsx file
write_xlsx(startClones, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "ProportionOfDay0ClonesThatAreCrossReactive.xlsx"))



#####
#Make a formal plot for clonal diversity using the Simpson index
# placeholder <- metadata[metadata$ProtoOmi != "Proto-Omi+",] %>%
#   group_by(Booster, adj.ProtoOmi, Timepoint, clone_subject_id) %>%
#   summarize(n= n())
# 
# groups <- unique(placeholder$Booster)
# df <- data.frame()
# 
# for(j in groups){
#   place <- placeholder[placeholder$Booster == j,]
#   cloneList <- split(place, list(place$Booster, place$Timepoint, place$adj.ProtoOmi))
#   
#   hold <- list()
#   for(i in 1:length(cloneList)){
#     hold[[i]] <- as.vector(cloneList[[i]]$n)
#   }
#   
#   names(hold) <- names(cloneList)
#   
#   diversity <- iNEXT(hold, q=2, datatype="abundance", endpoint=1500)
#   
#   df <- bind_rows(df, ggplot2::fortify(diversity, type=1))
#   
# }
# 
# df$Booster <- str_extract(df$Assemblage, "[A-Za-z ]+(?=\\.)")
# df$Timepoint <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9 ]+(?=\\.)")
# df$adj.ProtoOmi <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9+-]+\\Z")
# df$Timepoint <- factor(df$Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180"))
# 
# 
# df.point <- df[which(df$Method=="Observed"),]
# df.line <- df[which(df$Method!="Observed"),]
# df.line$Method <- factor(df.line$Method, c("Rarefaction", "Extrapolation"), c("Rarefaction", "Extrapolation"))
# 
# pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalDiversity_SimpsonIndex_Uninfected.pdf"),width=12, height = 12)
# ggplot(df, aes(x=x, y=y, color = Booster))+
#   geom_point(aes(shape = Booster), size=5, data=df.point)+
#   geom_line(aes(linetype=Method), lwd=1.5, data=df.line)+
#   geom_ribbon(aes(ymin=y.lwr, ymax=y.upr, fill= Booster, color=NULL), alpha=0.2)+
#   facet_grid(cols=vars(Timepoint), rows=vars(adj.ProtoOmi))+
#   labs(x="Sample Size of B Cell Population", y="Clonal Diversity")+
#   theme_classic()+
#   theme(legend.position = "bottom",
#         legend.title=element_blank(),
#         text=element_text(size=18),
#         legend.box="vertical",
#         axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# dev.off()
#####

#############
#############
##Statistics
#############
#############
#####
#Testing differences in VH usage between the three groups
stats <- metadata %>%
          group_by(Booster, Subject, v_call) %>%
          summarise(n= n()) %>%
          mutate(Proportion = n / sum(n)) %>%
          select(Booster, Subject, v_call, Proportion) %>%
          pivot_wider(names_from = v_call, values_from = Proportion)

stats[is.na(stats)] <- 0

#write a csv so that Flavio can do the analysis on Prism
write.csv(stats, here::here("04_Analysis", "data_objects","06_repertoire_analysis", "VHUsageProportions_WideFormat_PerIndividual.csv"))

#####
#Testing SHM between day 0 and 180 and by boost group
stats <- metadata[metadata$Timepoint %in% c("Day 0", "Day 180"),] %>%
  group_by(clone_subject_id, Booster, Timepoint) %>%
  summarize(meanSHM = mean(mu_freq),
            n = n())








######
#Construction yard
######
##########
#Here I am writing a script to color by both specificity and day first observed
# metadata$SubjGroup <- paste0(metadata$Booster, ", ", metadata$Subject)
# df <- metadata[metadata$adj.ProtoOmi != "Proto-Omi+",]
# 
# pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "donutplot_test_DarkColorsCrossReactive_GreenDay0_BlueCrossTimepoint.pdf"))
# for(i in unique(df$SubjGroup)){
#   placeholder <- df[df$SubjGroup== i,]
#   
#   nonSinglets <- unique(placeholder$clone_subject_id[duplicated(placeholder$clone_subject_id) | duplicated(placeholder$clone_subject_id, fromLast=T)])
#   placeholder$CloneStatus <- ifelse(placeholder$clone_subject_id %in% nonSinglets, placeholder$clone_subject_id, "Singlet")          
#   
#   placeholder$Timepoint <- ifelse(placeholder$Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", placeholder$Timepoint)
#   
#   calcs <- placeholder %>%
#     group_by(CloneStatus, Timepoint) %>%
#     summarize(n = n()) %>%
#     mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
#                            length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Day 0 Expanded",
#                            length(unique(Timepoint)) > 1 ~ "Expanded",
#                            TRUE ~ "Single Timepoint"))
#   
#   placeholder <- placeholder %>%
#     group_by(Timepoint, CloneStatus) %>%
#     summarize(n= n()) %>%
#     mutate(Proportion = n / sum(n),
#            Total = sum(n),
#            CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
#            adj.CloneStatus = case_when( CloneStatus == "Singlet" ~ "Singlet",
#                                         CloneStatus %in% calcs$CloneStatus[calcs$lab == "Day 0 Expanded"] ~ "Day 0 Expanded",
#                                         CloneStatus %in% calcs$CloneStatus[calcs$lab == "Expanded"] ~ "Expanded",
#                                         CloneStatus %in% calcs$CloneStatus[calcs$lab == "Single Timepoint"] ~ "Single Timepoint"),
#            adj.CloneStatus = fct(adj.CloneStatus, levels = c("Day 0 Expanded", "Expanded", "Single Timepoint", "Singlet")),
#            Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90/180")),
#            adj.ProtoOmi = df$adj.ProtoOmi[match(CloneStatus, df$clone_subject_id)],
#            adj.ProtoOmi = ifelse(is.na(adj.ProtoOmi), "", adj.ProtoOmi),
#            adj.CloneStatusSpec = ifelse(adj.CloneStatus == "Singlet", "Singlet", paste0(adj.CloneStatus," ", adj.ProtoOmi)))%>%
#     arrange(adj.CloneStatusSpec)%>%
#     mutate(ymax = cumsum(Proportion),
#            ymin = c(0, head(ymax, n=-1)))
#   
#   label <- c()
#   label[1] <- unique(placeholder$Total[placeholder$Timepoint == "Day 0"])
#   label[2] <- unique(placeholder$Total[placeholder$Timepoint == "Day 15"])
#   label[3] <- unique(placeholder$Total[placeholder$Timepoint == "Day 90/180"])
#   time <- c("Day 0", "Day 15", "Day 90/180")
#   dat_text <- data.frame(label = label, Timepoint = time)
#   dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90/180"))
#   
#   p <- ggplot(placeholder)+
#     geom_rect(color= "black", linewidth=0.2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=adj.CloneStatusSpec))+
#     coord_polar(theta="y")+
#     xlim(c(2,4))+
#     scale_fill_manual(values = c("Singlet" = "#FFFFFF",
#                                  "Day 0 Expanded Proto+Omi+" = "green4", "Day 0 Expanded Proto+Omi-" = "green",
#                                  "Expanded Proto+Omi+" = "navy", "Expanded Proto+Omi-" = "blue",
#                                  "Single Timepoint Proto+Omi+" = "gray25", "Single Timepoint Proto+Omi-" = "gray"))+
#     ggtitle(paste0("Subject, Group: ",i))+
#     guides(fill = "none")+
#     facet_grid(cols=vars(Timepoint))+
#     theme_void()+
#     theme(plot.title = element_text(hjust=0.5))+
#     geom_text(data = dat_text,
#               mapping = aes(x=-Inf, y=-Inf, label = label),
#               hjust = 0.5,
#               vjust = 0.5)
#   
#   print(p)
# }
# dev.off()
# rm(p)
# }
# rm(splitDF)
# rm(p)
#####

#####
#this code calculates and projects divesity of response per booster per specificity - consider adding a section that specifically looks at
#how the two compare at each timepoint side by side
# placeholder <- metadata[metadata$ProtoOmi != "Proto-Omi+",] %>%
#   group_by(Booster, adj.ProtoOmi, Timepoint, clone_subject_id) %>%
#   summarize(n= n())
# 
# groups <- unique(placeholder$Booster)
# df <- data.frame()
# 
# for(j in groups){
#   place <- placeholder[placeholder$Booster == j,]
#   cloneList <- split(place, list(place$Booster, place$Timepoint, place$adj.ProtoOmi))
#   
#   hold <- list()
#   for(i in 1:length(cloneList)){
#     hold[[i]] <- as.vector(cloneList[[i]]$n)
#   }
#   
#   names(hold) <- names(cloneList)
#   
#   diversity <- iNEXT(hold, q=2, datatype="abundance", endpoint=1500)
#   
#   #pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", paste0("ClonalDiversity_",j,".pdf")),width=25, height = 6)
#   #ggiNEXT(diversity, type=1, facet.var = "Assemblage", color.var = "Assemblage") + theme_classic() + theme(axis.text.x = element_text(angle=90, hjust = 1, vjust=0.5))
#   #dev.off()
#   
#   df <- bind_rows(df, ggplot2::fortify(diversity, type=1))
#   
# }
# 
# df$Booster <- str_extract(df$Assemblage, "[A-Za-z ]+(?=\\.)")
# df$Timepoint <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9 ]+(?=\\.)")
# df$adj.ProtoOmi <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9+-]+\\Z")
# 
# df.point <- df[which(df$Method=="Observed"),]
# df.line <- df[which()]
#####

#####
#For infected
# metadataInf <- seuObj@meta.data
# metadataInf <- metadataInf[!metadataInf$adj.ProtoOmi %in% c("Proto-Omi+", "Proto-Omi-"),]
# 
# placeholder <- metadataInf %>%
#   group_by(Booster, Infection, adj.ProtoOmi, Timepoint, clone_subject_id) %>%
#   summarize(n= n()) %>%
#   mutate(BoosterInf = paste0(Booster,"_",Infection)) %>%
#   ungroup()
# 
# groups <- unique(placeholder$BoosterInf)
# df <- data.frame()
# 
# for(j in groups){
#   place <- placeholder[placeholder$BoosterInf == j,]
#   place <- droplevels(place)
#   cloneList <- split(place, list(place$BoosterInf, place$Timepoint, place$adj.ProtoOmi))
#   
#   hold <- list()
#   for(i in 1:length(cloneList)){
#     hold[[i]] <- as.vector(cloneList[[i]]$n)
#   }
#   
#   names(hold) <- names(cloneList)
#   
#   diversity <- iNEXT(hold, q=2, datatype="abundance", endpoint=1500)
#   
#   df <- bind_rows(df, ggplot2::fortify(diversity, type=1))
#   
# }
# 
# df$BoosterInf <- str_extract(df$Assemblage, "[A-Za-z _]+(?=\\.)")
# df$Timepoint <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9 ]+(?=\\.)")
# df$adj.ProtoOmi <- str_extract(df$Assemblage, "(?<=\\.)[A-Za-z0-9+-]+\\Z")
# df$Timepoint <- factor(df$Timepoint, levels=c("Day 0", "Day 15", "Day 90", "Day 180"))
# 
# df.point <- df[which(df$Method=="Observed"),]
# df.line <- df[which(df$Method!="Observed"),]
# df.line$Method <- factor(df.line$Method, c("Rarefaction", "Extrapolation"), c("Rarefaction", "Extrapolation"))
# 
# pdf(file = here::here("04_Analysis", "plots", "06_repertoire_analysis", "ClonalDiversity_SimpsonIndex_Infected.pdf"),width=12, height = 12)
# ggplot(df, aes(x=x, y=y, color = adj.ProtoOmi))+
#   geom_point(aes(shape = adj.ProtoOmi), size=5, data=df.point)+
#   geom_line(aes(linetype=Method), lwd=1.5, data=df.line)+
#   geom_ribbon(aes(ymin=y.lwr, ymax=y.upr, fill= adj.ProtoOmi, color=NULL), alpha=0.2)+
#   facet_grid(cols=vars(Timepoint), rows=vars(BoosterInf))+
#   labs(x="Sample Size of B Cell Population", y="Clonal Diversity")+
#   theme_classic()+
#   theme(legend.position = "bottom",
#         legend.title=element_blank(),
#         text=element_text(size=18),
#         legend.box="vertical",
#         axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
# dev.off()

#####