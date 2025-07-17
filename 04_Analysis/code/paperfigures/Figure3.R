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
library(ggalluvial)

set.seed(1)

#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#2AB673", 
               "Prototype + Omicron BA.1 mRNA" = "#1D75BC",
               "Prototype mRNA" = "#FBB042")

immunogenColors <- c("Prototype" = "#FBB042",
                     "BA.1 And Prototype" = "#1D75BC",
                     "BA.1" = "#2AB673")

shortColors <- c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 #"Acute Activated" = "#F46D43",
                 "Acute Activated" = "#f08665",
                 "Intermediate" = "#E6F598",
                 "Resting IgG" = "limegreen",
                 "Resting IgA" = "#3288BD",
                 "Plasmablast-like" = "#6f2da8",
                 "Naive" = "white")

#load in the data
seuObjNaive <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObjNaive %>% filter(ClusterLabel != "Naive")

df <- seuObj@meta.data
df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"),
         Booster = str_replace(Booster, "Omicron", "BA.1"),
         OfficialBooster = factor(OfficialBooster, levels = c("Prototype mRNA", "Prototype + Omicron BA.1 mRNA", "Omicron BA.1 mRNA")),
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")))

naiveDF <- seuObjNaive@meta.data
embeddings <- as.data.frame(Embeddings(seuObjNaive,reduction = "harmony.wnn.umap"))
naiveDF$UMAP_1 <- embeddings$harmonywnnUMAP_1[match(rownames(embeddings), rownames(naiveDF))]
naiveDF$UMAP_2 <- embeddings$harmonywnnUMAP_2[match(rownames(embeddings), rownames(naiveDF))]

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
#####

#####
#create a table for intro graphic
table(naiveDF$Booster, naiveDF$Timepoint)
#####

######
#Figure 3: showing clustering results of all cells
ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = ClusterLabel))+
  geom_point(size = 0.9, shape=21, stroke = 0.1)+
  scale_fill_manual(values = shortColors)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  theme(legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(0.1, "line"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3a_UMAP.png"),width = 2.9, height = 2.7, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3a_UMAP.svg"),width = 2.9, height = 2.7, units = "in")
dev.off()
#####

#####
#Make an alluvial plot showing phenotypes over time
stats <- naiveDF %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
group_by(Timepoint, ClusterLabel) %>%
  mutate(ClusterLabel = factor(ClusterLabel, levels = c("Acute Activated", "Intermediate", "Resting IgG", "Resting IgA","Atypical", "Plasmablast-like", "Naive"))) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum(linewidth = 0.4)+
  scale_fill_manual(values = shortColors)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  theme(legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "lines"),
        axis.title.y = element_text(size=7.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7.5,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(size=6, face = "bold"),
        panel.spacing = unit(0.4, "lines"))+
  guides(shape = guide_legend(override.aes = list(size=0.4)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "GeneralAlluvialPlot.png"),width = 3.2, height = 2.2, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "GeneralAlluvialPlot.svg"),width = 3.2, height = 2.2, units = "in", device = "svg")
dev.off()
#####

#####
#Make a heatmap and violin plot of phenotypes
phenoColors <- c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 "Acute Activated" = "#f08665",
                 "Intermediate" = "#E6F598",
                 "Resting IgG" = "limegreen",
                 "Resting IgA" = "#3288BD",
                 "Plasmablast-like" = "#6f2da8",
                 "Naive" = "black") #you can't scale transparency when a color is white LMAO

#Heatmap showing IGC usage 
stats <- naiveDF %>% filter(!is.na(c_call)) %>%
  group_by(ClusterLabel, c_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(ClusterLabel, c_call, fill = list(Proportion = 0))

#add in cluster label without memory- change color scheme to short color
ggplot(stats, aes(x=c_call, y=ClusterLabel, alpha=Proportion))+
  #geom_tile(color="white", linewidth=1, fill="dodgerblue4")+
  geom_tile(color="white", linewidth=1, aes(fill=ClusterLabel))+
  scale_alpha(range = c(0.1, 1.4), limits = c(0, 0.8))+
  scale_fill_manual(values = phenoColors, guide = "none")+
  coord_equal()+
  #scale_y_discrete(limits=rev(c("Week 0","Week 1","Week 4", "Week 16", "Week 17", "Week 18", "Week 20", "Week 28", "Week 40")))+
  scale_x_discrete(limits = c("IGHM", "IGHD", "IGHG1", "IGHG3", "IGHG2", "IGHG4", "IGHA1", "IGHA2"), position = "top")+
  scale_y_discrete(limits = rev(c("Atypical", "Acute Activated", "Intermediate","Resting IgG", "Resting IgA", "Plasmablast-like", "Naive")), position = "left")+
  #labs(y="Cluster")+
  #scale_fill_brewer(palette="Blues")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 270, vjust = 0, hjust=0.5, size=8), 
        axis.text.y = element_text(size=6),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        #legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        legend.position = "left",
        legend.key.size = unit(0.1, "lines"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_FCProportionalUsage.png"),width = 2.7, height = 3.3, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_FCProportionalUsage.svg"),width = 2.7, height = 3.3, units = "in")
dev.off()

##alternative since stacking violins is difficult
#important protein markers to show:
#IgA, IgG, CD21, CD27, CD11c, CD62L, FCRL5
#Other interesting ones:
#CD72, CD95, CD45-RB (see https://www.cell.com/immunity/pdf/S1074-7613(20)30268-5.pdf),
protein <- as.data.frame(t(seuObjNaive@assays[["Prot"]]@scale.data))
dfExtended <- cbind(naiveDF, protein) %>% select(CELL, ClusterLabel, colnames(protein)) %>%
  pivot_longer(!c(CELL, ClusterLabel), names_to = "marker", values_to="expression") %>%
  mutate(marker = gsub("P-", "", marker)) %>%
  #filter(marker %in% c("IgA", "IgG", "IgM", "CD21", "CD27", "CD11c", "CD62L", "CD85j", "FCRL5", "CD72", "CD95", "CD45RB")) %>%
  #mutate(marker = factor(marker, levels = c("IgG", "IgA", "IgM", "CD21", "CD27", "CD62L", "CD85j", "CD11c", "FCRL5", "CD72", "CD95", "CD45RB")))
  filter(marker %in% c("IgA", "IgG", "IgM", "IgD", "CD21", "CD71" ,"CD72", "CD27", "CD11c", "FCRL5","CD45RB")) %>%
  mutate(marker = factor(marker, levels = c("IgD", "IgM", "IgA", "IgG", "CD21", "CD27", "CD71","CD72", "CD11c", "FCRL5","CD45RB")))

#markers to include:
#CD21 bc she's a baddie
#CD11c bc it's a classic atypical marker
#CD45RB shown to be downregulated on atypicals in the context of SARS-CoV2
ggplot(dfExtended, aes(y= ClusterLabel, x = expression))+
  geom_violin(aes(fill = ClusterLabel), linewidth = 0.2)+
  scale_fill_manual(values = phenoColors)+
  scale_x_continuous(position = "bottom", n.breaks = 4)+
  scale_y_discrete(limits = rev(c("Atypical","Acute Activated", "Intermediate","Resting IgG", "Resting IgA", "Plasmablast-like", "Naive")), position = "left")+
  facet_grid(cols = vars(marker), scales= "free_x")+
  annotate("rect", xmin=-Inf, xmax = Inf, ymin= -Inf, ymax = Inf, fill = NA, color = "black")+
  geom_vline(xintercept = -Inf, linewidth = 0.8)+
  theme_classic()+
  theme(legend.text = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        #axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=17),
        strip.text.x = element_text(angle = 270, size = 9, vjust = 0.5, hjust = 1),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 7, angle = 270, vjust = 0.5),
        panel.spacing =unit(0.2, "lines"),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_ProteinViolinPlot_reduced.png"),width = 4.2, height = 2.2, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3_ProteinViolinPlot_reduced.svg"),width = 4.2, height = 2.2, units = "in")
dev.off()
#####

#####
#Make a barplot showing general CITESeq specificities
######
stats <- df %>%
  group_by(adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         adj.ProtoOmi = factor(adj.ProtoOmi, levels = c("Proto+Omi+",
                                                        "Proto+Omi-","Proto-Omi+")))

ggplot(stats, aes(x = 1, y = Proportion))+
  geom_bar(stat = "identity", position = "stack", aes(fill = adj.ProtoOmi))+
  scale_fill_manual(values = c("Proto+Omi+" = "black",
                               "Proto+Omi-" = "gray50",
                               "Proto-Omi+" = "gray80"))+
  ylab("Proportion")+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "TotalCITESeqBarplot_Specificities.png"), width = 2.1, height = 2)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "TotalCITESeqBarplot_Specificities.svg"), width = 2.1, height = 2)
######

#####
#Show correlation between citeseq and flow data
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

#########
#Make an alluvial plot and split by specificity
stats <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), ClusterLabel != "Plasmablast-like", adj.ProtoOmi != "Proto-Omi+") %>%
  mutate(adj.ProtoOmi = case_when(adj.ProtoOmi == "Proto+Omi+" ~ "Cross-Reactive",
                                  adj.ProtoOmi == "Proto+Omi-" ~ "Prototype-Specific")) %>%
  group_by(Booster, adj.ProtoOmi ,Timepoint, ClusterLabel) %>%
  mutate(ClusterLabel = factor(ClusterLabel, levels = c("Acute Activated", "Intermediate", "Resting IgG", "Resting IgA","Atypical", "Plasmablast-like", "Naive"))) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum(linewidth = 0.4)+
  scale_fill_manual(values = shortColors)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), labeller = label_wrap_gen(width = 15))+
  theme(legend.text = element_text(size=7),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(0.4, "lines"),
        axis.title.y = element_text(size=7.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7.5,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7),
        strip.background = element_blank(),
        strip.text = element_text(size=6, face = "bold"),
        panel.spacing = unit(0.4, "lines"))+
  guides(shape = guide_legend(override.aes = list(size=0.4)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3b_AlluvialPlotOfPopulations_splitbyspecificity.png"),width = 4.2, height = 2.6, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "Figure3b_AlluvialPlotOfPopulations_splitbyspecificity.svg"),width = 4.2, height = 2.6, units = "in", device = "svg")
dev.off()
########

######
#Compare activation by specificity
########
calcs <- df %>% filter(Timepoint %in% c("Day 0", "Day 15"), adj.ProtoOmi != "Proto-Omi+") %>%
  mutate(ClusterLabel = case_when(ClusterLabel %in% c("Intermediate", "Acute Activated") ~ "Activated",
                                  TRUE ~ ClusterLabel)) %>%
  mutate(ClusterSpec = paste0(ClusterLabel,"_", adj.ProtoOmi)) %>%
  group_by(Booster, Subject, Timepoint, ClusterSpec) %>%
  summarize(n = n()) %>%
  group_by(Booster)%>%
  complete(Subject, Timepoint, ClusterSpec, fill = list(n = 0)) %>%
  group_by(Booster, Subject, Timepoint)%>%
  mutate(Prop = n / sum(n)) %>%
  mutate(ClusterLabel = str_extract(ClusterSpec, "Activated"),
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Activated Cross-Reactive", "Activated Prototype-Specific")) %>%
  filter(!is.na(ClusterLabel))

#plot out the differences per group
#prototype
calcs %>% filter(Booster == "Prototype") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 0.9)+
  ylab("Proportion of Total Cells")+
  ggtitle("Prototype Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(adj.ProtoOmi))+
  scale_fill_manual(values = c("Activated Prototype-Specific" = "lightgray",
                               "Activated Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype-Specific" = "lightgray",
                                "Activated Cross-Reactive" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_PrototypeBoost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_PrototypeBoost.svg"), width = 1.6, height = 1.8)

#proto/BA1
calcs %>% filter(Booster == "BA.1 And Prototype") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 0.9)+
  ylab("Proportion of Total Cells")+
  ylim(0,1)+
  ggtitle("Prototype/BA.1 Boost")+
  facet_grid(cols = vars(adj.ProtoOmi))+
  scale_fill_manual(values = c("Activated Prototype-Specific" = "lightgray",
                               "Activated Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype-Specific" = "lightgray",
                                "Activated Cross-Reactive" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_BivalentBoost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_BivalentBoost.svg"), width = 1.6, height = 1.8)

#BA.1 vaccine
calcs %>% filter(Booster == "BA.1") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 0.9)+
  ylab("Proportion of Total Cells")+
  ggtitle("BA.1 Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(adj.ProtoOmi))+
  scale_fill_manual(values = c("Activated Prototype-Specific" = "lightgray",
                               "Activated Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype-Specific" = "lightgray",
                                "Activated Cross-Reactive" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_BA1Boost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 3", "ProportionTotalCells_Activated_BA1Boost.svg"), width = 1.6, height = 1.8)
######