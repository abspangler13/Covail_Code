library(ggplot2)
library(dplyr)
library(here)
library(Seurat)
library(readxl)
library(tidyseurat)
library(stringr)
library(alakazam)
library(scoper)
library(writexl)
library(stats)
library(ggseqlogo)
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

######
#make VH gene plot using significant specificity genes
vh <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
  group_by(Subject, adj.ProtoOmi, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup()%>%
  complete(Subject, adj.ProtoOmi, v_call, fill = list(n = 0, Proportion = 0)) %>%
  group_by(adj.ProtoOmi, v_call) %>%
  summarize(median = mean(Proportion)) %>%
  mutate(adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Omicron Cross-Reactive",
                                  adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Specific"))

#plot now
ggplot(vh, aes(x = v_call, y = adj.ProtoOmi, alpha = median))+
  geom_tile(color = "white", fill = "#002D04", linewidth = 1)+
  #scale_fill_manual(values = c("Omicron Cross-Reactive" = "#1AAFBC",
  #                             "Prototype-Specific" = "#80634C"))+
  scale_alpha(limits= c(0,0.16),range = c(0.01, 1))+
  ylab("Specificity")+
  xlab("VH Gene")+
  scale_y_discrete(limits = c("Omicron Cross-Reactive", "Prototype-Specific"), labels = c("Omicron\nCross-Reactive", "Prototype-Specific"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        legend.key.size = unit(0.3, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "VHHeatmap_CrossvsProto.png"), width = 7.4, height = 0.9, dpi = 1200)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "VHHeatmap_CrossvsProto.svg"), width = 7.4, height = 0.9)

###test statistically?
vh <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
  group_by(Subject, adj.ProtoOmi, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(Subject, adj.ProtoOmi, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "adj.ProtoOmi", values_from = "Proportion")

#do stats and compare VH genes
uniqueVH <- unique(vh$v_call)

pvals <- c()
vhGene <- c()
for(i in uniqueVH){
  vhGene <- append(vhGene, i)
  vhFiltered <- vh %>% filter(v_call == i)
  
  pvals <- append(pvals, wilcox.test(vhFiltered$`Proto+Omi+`, vhFiltered$`Proto+Omi-`, paired= TRUE)$p.value)
}

pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))


#plotting out specific VH genes?
gene <- unique(pvDF$vhGene[pvDF$pvals < 0.05]) #choose a gene
vhSingle <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
  group_by(Subject, adj.ProtoOmi, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(Subject, adj.ProtoOmi, v_call, fill = list(Proportion = 0)) %>%
  filter(v_call %in% gene) %>%
  mutate(adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Cross-Reactive",
                                  adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Only"))

ggplot(vhSingle, aes(x = adj.ProtoOmi, y = Proportion, fill = adj.ProtoOmi))+
  #geom_point(shape = 21, position = position_jitter(width = 0.3))+
  geom_line(aes(group = Subject), alpha = 0.5, linewidth = 0.1)+
  geom_point(shape = 21, size = 0.6)+
  facet_grid(cols = vars(v_call))+
  scale_fill_manual(values = c("Cross-Reactive" = "#ADEFD1",
                               "Prototype-Only" = "#00203F"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "SignificantVHGenes.png"), width = 4, height = 1.5)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "SignificantVHGenes.svg"), width = 4, height = 1.5)
######

######
#identify and plot out public clones and their specificities
vh <- df %>% filter(adj.ProtoOmi != "Proto-Omi+") %>%
  group_by(Subject, adj.ProtoOmi, v_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  select(!n) %>%
  ungroup() %>%
  complete(Subject, adj.ProtoOmi, v_call, fill = list(Proportion = 0)) %>%
  pivot_wider(names_from = "adj.ProtoOmi", values_from = "Proportion")

#do stats and compare VH genes
uniqueVH <- unique(vh$v_call)

pvals <- c()
vhGene <- c()
for(i in uniqueVH){
  vhGene <- append(vhGene, i)
  vhFiltered <- vh %>% filter(v_call == i)
  
  pvals <- append(pvals, wilcox.test(vhFiltered$`Proto+Omi+`, vhFiltered$`Proto+Omi-`, paired= TRUE)$p.value)
}

pvDF <- data.frame(vhGene, pvals) %>% mutate(pvals = p.adjust(pvals, method = "bonferroni"))

#slim down the data to just the DE VH genes
slimmedDF <- df %>% filter(v_call %in% pvDF$vhGene[pvDF$pvals <= 0.05], adj.ProtoOmi != "Proto-Omi+")%>%
  select(!clone_id)

set.seed(1)

results <- hierarchicalClones(slimmedDF,
                              method = "nt",
                              threshold = 0.20,
                              only_heavy = TRUE, split_light = TRUE,
                              summarize_clones = FALSE)

summary <- results %>% group_by(clone_id) %>% mutate(Public = length(unique(Subject))) %>% filter(Public >= 10)

#check <- summary %>% group_by(clone_id) %>% filter(n() > 20) %>% select(clone_id, junction_aa)

###########first make the donut plots showing the proportion of specificities per clone id
summaryDonut <- summary %>% group_by(clone_id, v_call, adj.ProtoOmi) %>% summarize(n = n(), Public = mean(Public)) %>% 
  mutate(sum = sum(n), Proportion = n / sum, ymax = cumsum(Proportion), ymin =c(0, head(ymax, n=-1)), 
         adj.ProtoOmi = ifelse(adj.ProtoOmi == "Proto+Omi+", "Cross-Reactive", "Prototype-Only")) %>%
  filter(sum > 25, clone_id %in% c(2433, 3480, 3111, 1817, 1893)) %>%
  mutate(clone_id = paste0("Public", clone_id, ", ", v_call),
         clone_id = factor(clone_id, levels = c("Public2433, IGHV2-70", 
                                                "Public3480, IGHV3-30", "Public3111, IGHV3-30", 
                                                "Public1817, IGHV5-51", "Public1893, IGHV5-51")))

ggplot(summaryDonut)+
  geom_rect(color= "black", linewidth=0.1, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=adj.ProtoOmi))+
  coord_polar(theta="y")+
  xlim(c(2,4))+
  scale_fill_manual(values = c("Cross-Reactive" = "#ADEFD1",
                               "Prototype-Only" = "#00203F"))+
  geom_text(#data = dat_text,
    mapping = aes(x=-Inf, y=-Inf, label = sum),
    hjust = 0.5,
    vjust = 0.5,
    size = 3)+
  facet_grid(rows=vars(clone_id), labeller = label_wrap_gen(10), switch="y")+
  guides(fill = guide_legend(ncol= 1))+
  theme_void()+
  theme(plot.title = element_text(size = 3, hjust = 0.5, vjust = 0.5),
        strip.text.x = element_text(size = 6),
        #panel.spacing = unit(-0.49999, "lines"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        panel.spacing.x = unit(1, "lines"),
        strip.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "PublicClones_SpecificityLabels.png"), width = 2.3, height = 4.3)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "PublicClones_SpecificityLabels.svg"), width = 2.3, height = 4.3)

#######make logo plots for each clonal group
seqL <- summary %>% group_by(clone_id) %>%  filter(clone_id %in% c(2433, 3480, 3111, 1817, 1893)) %>% select(clone_id, junction_aa) %>% ungroup()
cloneList <- c("2433", "3480", "3111", "1817", "1893")
seqList <- list()
for(i in cloneList) {
  seqList[[i]] <- seqL$junction_aa[seqL$clone_id == i]
}

plotList <- list()
for(j in names(seqList)){
  p <- ggseqlogo(seqList[j], seq_type = "aa", method = "prob")+
    theme(text = element_text(size = 6),
          legend.position = "none",
          plot.margin = margin(0, 0, 0, 0))
  
  plotList[[j]] <- p
}

g <- arrangeGrob(grobs = plotList, nrow = length(plotList), ncol=1)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "PublicClone_Logoplots.png"), width = 4, height = 4, dpi = 1500)
ggsave(g, filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "PublicClone_Logoplots.svg"), width = 4, height = 4, dpi = 1500)
######

######
#calculate diversity
diversityDF <- df %>% filter(Infection == "N") %>%
  #filter(Timepoint %in% c("Day 0", "Day 15"), adj.ProtoOmi != "Proto-Omi+") %>%
  filter(adj.ProtoOmi != "Proto-Omi+") %>%
  mutate(BoosterDay = paste(OfficialBooster, Subject, Timepoint, adj.ProtoOmi))

###method 1: calculating using the alphadiversity function
curve <- alakazam::alphaDiversity(diversityDF, group="BoosterDay", min_q=2, max_q=2, min_n=20) #calculate bootstrapped diversity estimate at q=2 (Simpson)

calculatedDiversity <- curve@diversity %>% filter(q == 2) %>% #convert data to graphable form
  mutate(Booster = str_remove(BoosterDay, " [0-9]+ Day [0-9]+ Proto[+-]Omi[+-]"),
         Timepoint = str_extract(BoosterDay, "Day [0-9]+"),
         adj.ProtoOmi = str_extract(BoosterDay, "Proto[+-]Omi[+-]"),
         Donor = str_extract(BoosterDay, "(?<=mRNA )[0-9]+")
  ) %>%
  group_by(Donor, adj.ProtoOmi) %>%
  mutate(AllTimes = length(unique(Timepoint)) < 2,
         adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Omicron Cross-Reactive",
                                  adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Specific")) %>%
  filter(!AllTimes)

ggplot(calculatedDiversity, aes(x = Timepoint, y = d))+
  geom_point(aes(fill = Booster), shape =21, size = 1)+
  geom_line(aes(group = Donor, color = Booster), linewidth = 0.3, alpha = 0.5)+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+  
  scale_x_discrete(limits= c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(rows = vars(Booster), cols = vars(adj.ProtoOmi), axes = "all")+
  ylim(8, 21)+
  ylab("Simpson's Index")+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "diversity_diversityCurve_perdonor.png"),width = 2.4, height = 3, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "diversity_diversityCurve_perdonor.svg"),width = 2.4, height = 3, units = "in", dpi = 1200)
dev.off()

#write a table
calculatedDiversity %>% select(!c(d_sd, d_lower, d_upper, e, e_lower, e_upper)) %>% select(Booster, Donor, Timepoint, adj.ProtoOmi, q, d) %>% pivot_wider(names_from =  Timepoint, values_from = d) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S4", "diversity_perindividual_bootstrap.xlsx"))
######

######
#Look at % VH mutation for prototype-specific cells
shm <- df %>% filter(adj.ProtoOmi == "Proto+Omi-", Infection == "N") %>% mutate(mu_freq = mu_freq * 100)

ggplot(shm, aes(x = Timepoint, y = mu_freq, fill = OfficialBooster))+
  geom_violin()+
  ylab("% VH Mutation")+
  geom_boxplot(width=0.2, outlier.size = 0)+
  facet_grid(cols = vars(OfficialBooster), labeller = label_wrap_gen(width = 15))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors)+
  ggtitle("Prototype+ Omicron-")+
  theme_classic()+
  ylim(0, 12.5)+
  theme(title = element_text(size = 8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 6),
        strip.background = element_blank(),
        strip.text = element_text(size = 7, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "Figure3_SHM_protospecific.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "Figure3_SHM_protospecific.svg"),width = 3, height = 1.8, units = "in")
dev.off()

#do stats
shm %>% select(Timepoint, Booster, mu_freq) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S4", "mu_freq_OverTime_prototypeonly.xlsx"))
#####

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
  ggtitle("Prototype+ Omicron-")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols= vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 7),
        strip.text = element_text(size = 8),
        legend.position = "none",
        plot.title = element_text(size = 6, hjust = 0.5),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "SHMOverall_by_clone_proto.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S4", "SHMOverall_by_clone_proto.svg"),width = 3, height = 1.8)
dev.off()

#write a stats sheet
shm %>%
  group_by(Booster, clone_subject_id, Timepoint) %>%
  summarize(n = n(),
            mean = mean(mu_freq)) %>%
  group_by(clone_subject_id) %>%
  mutate(Present = length(unique(Timepoint)) == 4) %>%
  filter(Present) %>%
  group_by(Booster,clone_subject_id)%>% select(!c(n, Present)) %>%
  pivot_wider(names_from = Timepoint, values_from = mean) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S4", "SHM_perClonalLineage_proto.xlsx"))
#####

######
#Plot pre- vs post-infection SHM
crossDF <- df %>% filter(Infection == "Y") %>% filter(adj.ProtoOmi == "Proto+Omi+")

nonSinglets <- unique(crossDF$clone_subject_id[duplicated(crossDF$clone_subject_id) | duplicated(crossDF$clone_subject_id, fromLast=T)])
crossDF$CloneStatus <- ifelse(crossDF$clone_subject_id %in% nonSinglets, crossDF$clone_subject_id, "Singlet")

calcs <- crossDF %>%
  group_by(CloneStatus, InfectionRange, Timepoint) %>%
  summarize(n = n()) %>%
  group_by(CloneStatus, InfectionRange) %>%
  mutate(relative = case_when(unique(InfectionRange) == "Between Days 15-90" & Timepoint %in% c("Day 90", "Day 180") ~ "Post-infection",
                              unique(InfectionRange) == "Between Days 90-180" & Timepoint %in% c("Day 180") ~ "Post-infection",
                              TRUE ~ "Pre-infection")) %>%
  mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) & "Post-infection" %in% unique(relative) ~ "Expanded Pre-Vax, Post-Infection",
                         length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Expanded Pre-Vax",
                         length(unique(Timepoint)) > 1 & length(unique(relative)) == 2 ~ "Expanded Pre- and Post-Infection",
                         length(unique(Timepoint)) > 1 & length(unique(relative)) == 1 & "Pre-infection" %in% unique(relative) ~ "Expanded Pre-Infection",
                         length(unique(relative)) == 1 & "Post-infection" %in% unique(relative) ~ "Expanded Post-Infection",
                         length(unique(relative)) == 1 & !("Post-infection" %in% unique(relative)) ~ "Single Timepoint Pre-Infection",
                         TRUE ~ "CHECK"))

crossDF$CloneStatusRefined <- calcs$lab[match(crossDF$CloneStatus, calcs$CloneStatus)]
crossDF$CloneStatusRefined <- factor(crossDF$CloneStatusRefined, levels = c("Expanded Pre-Vax, Post-Infection",
                                                                            "Expanded Pre- and Post-Infection",
                                                                            "Expanded Post-Infection",
                                                                            "Expanded Pre-Vax",
                                                                            "Expanded Pre-Infection",
                                                                            "Single Timepoint Pre-Infection",
                                                                            "Singlet"))

compareDF <- crossDF %>% filter(CloneStatusRefined %in% c("Expanded Post-Infection", "Expanded Pre-Vax, Post-Infection",
                                                          "Expanded Pre-Infection", "Expanded Pre- and Post-Infection")) %>%
  mutate(CloneStatus2 = case_when(CloneStatusRefined == "Expanded Post-Infection" ~ "Post-Infection",
                                  TRUE ~ "Pre-Infection"))

#plot
ggplot(compareDF, aes(x = CloneStatus2, y = mu_freq*100, fill = CloneStatus2))+
  geom_violin()+
  geom_boxplot(outlier.size = 0, width = 0.1)+
  scale_x_discrete(limits = c("Pre-Infection", "Post-Infection"))+
  scale_fill_manual(values = c("Post-Infection" = "#FFAE42",
                               "Pre-Infection" = "#A8415B"))+
  ylab("% VH Mutation")+
  theme_classic()+
  theme(text = element_text(size = 9, color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust =1, vjust = 1, color = "black"),
        axis.line.x = element_line(linewidth = 0.2),
        axis.line.y = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2))
ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure S4", "Pre_Post_Infection_SHM.png"), width = 2, height =2.3)
ggsave(here::here("04_Analysis", "plots","paperfigures", "Figure S4", "Pre_Post_Infection_SHM.svg"), width = 2, height =2.3)
wilcox.test(compareDF$mu_freq[compareDF$CloneStatus2 == "Post-Infection"],
            compareDF$mu_freq[compareDF$CloneStatus2 == "Pre-Infection"], paired = FALSE)
######