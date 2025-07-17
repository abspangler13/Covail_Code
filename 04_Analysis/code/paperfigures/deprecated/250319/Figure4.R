library(Seurat)
library(dplyr)
library(tidyseurat)
library(ggplot2)
library(ggalluvial)
library(here)

#load in the data
seuObjNaive <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
seuObj <- seuObjNaive %>% filter(ClusterLabel != "Naive")
df <- seuObj@meta.data #creates a dataframe of the metadata

df <- df %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"),
         ClusterLabel = case_when(ClusterLabel == "Activated IgA Memory" ~ "Atypical IgA Memory",
                                  TRUE ~ ClusterLabel))


df$ClusterLabel <- factor(df$ClusterLabel, levels = c("AM1 (Activated)",
                                                      "AM2 (Intermediate)",
                                                      "AM3 (Atypical)",
                                                      "Atypical IgA Memory",
                                                      "Resting IgA Memory",
                                                      "Resting IgG Memory",
                                                      "Resting IgG Memory 2",
                                                      "Unclear"))

# df$ClusterLabel <- factor(df$ClusterLabel, levels = c("C1-Activated",
#                                                       "C2-Atypical",
#                                                       "C4-Activated IgA",
#                                                       "C3-Intermediate",
#                                                       "C6-Resting IgA",
#                                                       "C8-Resting IgG",
#                                                       "C9-Resting IgG 2",
#                                                       "Unclear"))

df <- df %>% mutate(ClusterLabelSep = ClusterLabel,
                    ClusterLabel = case_when(ClusterLabel == "Resting IgG Memory 2" ~ "Resting IgG Memory",
                                                     TRUE ~ ClusterLabel))

df$OfficialInfection <- ifelse(df$Infection == "Y", "Infected", "Uninfected")
df$OfficialInfection <- factor(df$OfficialInfection, levels = c("Uninfected",
                                                                "Infected"))

df$InfectionRange <- ifelse(is.na(df$InfectionRange), "Uninfected", df$InfectionRange)

#set colors
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#045275")

phenoColors <- c("AM1 (Activated)" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 "AM2 (Intermediate)" = "#F46D43",
                 "AM3 (Atypical)" = "#FDAE61",
                 "Atypical IgA Memory" = "#FEE08B",
                 "Resting IgA Memory" = "#E6F598",
                 "Resting IgG Memory" = "#ABDDA4",
                 #"Resting IgG Memory 2" = "#66C2A5",
                 "Unclear" = "#3288BD",
                 "Naive" = "#6f2da8")

# phenoColors <- c("C1-Activated" = "#D53E4F", #based on RColorBrewer Spectral Palette
#                  "C2-Atypical" = "#F46D43",
#                  "AM3 (Atypical)" = "#FDAE61",
#                  "Atypical IgA Memory" = "#FEE08B",
#                  "Resting IgA Memory" = "#E6F598",
#                  "Resting IgG Memory" = "#ABDDA4",
#                  "Resting IgG Memory 2" = "#66C2A5",
#                  "Unclear" = "#3288BD",
#                  "Naive" = "#6f2da8")



#####
#Figure 4: showing clustering results of all cells
#Let's work with data that includes naive cells just for now
naiveDF <- seuObjNaive@meta.data %>% mutate(ClusterLabelSep = ClusterLabel,
                                        ClusterLabel = case_when(ClusterLabel == "Resting IgG Memory 2" ~ "Resting IgG Memory",
                                                                 TRUE ~ ClusterLabel))

embeddings <- as.data.frame(Embeddings(seuObjNaive,reduction = "harmony.wnn.umap"))
naiveDF$ClusterLabel <- factor(naiveDF$ClusterLabel, levels = c("AM1 (Activated)",
                                                      "AM2 (Intermediate)",
                                                      "AM3 (Atypical)",
                                                      "Atypical IgA Memory",
                                                      "Resting IgA Memory",
                                                      "Resting IgG Memory",
                                                      "Resting IgG Memory 2",
                                                      "Unclear",
                                                      "Naive"))

naiveDF$UMAP_1 <- embeddings$harmonywnnUMAP_1[match(rownames(embeddings), rownames(naiveDF))]
naiveDF$UMAP_2 <- embeddings$harmonywnnUMAP_2[match(rownames(embeddings), rownames(naiveDF))]

#let's plot this out- first
ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = ClusterLabel))+
  geom_point(size = 1.1, shape=21, stroke = 0.15)+
  scale_fill_manual(values = phenoColors)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size=3)))+
  theme(legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_blank())
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4a_UMAP.png"),width = 3.5, height = 3, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 3b: showing phenotypic changes over time
stats <- df %>% filter(Infection == "N") %>%
  group_by(Booster, Timepoint, ClusterLabel) %>%
  summarize(n= n()) %>%
  mutate(Proportion = n / sum(n))

ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum(linewidth = 0.4)+
  scale_fill_brewer(palette = "Spectral")+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(width = 15))+
  theme(legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key.size = unit(0.6, "lines"),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7),
        #legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size=7))+
  guides(shape = guide_legend(override.aes = list(size=0.4)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4b_AlluvialPlotOfPopulations.png"),width = 4, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 4c: Plotting out specific AM compartments over time
#note from rory: not sure if I wanna actually plot this since AM2 is not actually AM2 from Sarah's H7 paper
stats <- df %>% filter(Infection == "N") %>%
          group_by(Timepoint, ClusterLabel) %>%
          summarize(n = n()) %>%
          mutate(Proportion = n / sum(n),
                 ShortenedLabel = substr(ClusterLabel, start = 1, stop = 3))

#let's plot out each 
#Each of the AMs
# ggplot(stats[stats$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)"),], aes(x = Timepoint, y = Proportion, fill = ShortenedLabel))+
#   geom_line(aes(group = ShortenedLabel, color = ShortenedLabel))+
#   geom_point(size = 2, shape = 21)+
#   ylab("Proportion")+
#   ylim(c(0,0.4))+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   scale_fill_manual(values = activatedColors)+
#   scale_color_manual(values = activatedColors)+
#   theme_classic()+
#   theme(legend.key.size = unit(0.1, 'cm'),
#     legend.text = element_text(size=7),
#     legend.title = element_text(size = 0),
#     axis.title.y = element_text(size=7),
#     axis.title.x = element_text(size=0),
#     axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#     axis.text.y = element_text(size = 7),
#     legend.margin=margin(0,0,0,0))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4c_AMPopulationsOverTime.png"),width = 2, height = 2, units = "in", device = "png", dpi = 600)
# dev.off()
#####

#####
#Fig 4
#Alluvial phenotype plots, splitting by Probe Specificity
# stats <- df %>% filter(Infection == "N" & adj.ProtoOmi != "Proto-Omi+") %>%
#   group_by(OfficialBooster, adj.ProtoOmi, Timepoint, ClusterLabel) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n))
# 
# #plot
# ggplot(stats, aes(y = Proportion, x= Timepoint, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
#   geom_flow()+
#   geom_stratum()+
#   scale_fill_brewer(palette = "Spectral")+
#   scale_y_continuous(expand = c(0,0))+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(OfficialBooster), rows = vars(adj.ProtoOmi), labeller=label_wrap_gen(width = 15))+
#   theme_classic()+
#   theme(legend.text = element_text(size = 6),
#         legend.key.size = unit(0.1, 'cm'),
#         legend.title = element_text(size = 7),
#         legend.margin=margin(0,0,0,0),
#         axis.title.y = element_text(size=8),
#         axis.title.x = element_text(size=0),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size = 8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 6, face = "bold"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_AMPopulationsOverTime_ByProbeSpec.png"),width = 5, height = 3.5, units = "in", device = "png", dpi = 600)
# dev.off()
#####

#####
#Probe specificity on the umap
# ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = adj.ProtoOmi))+
#   geom_point(size = 1.1, shape=21, stroke = 0.15)+
#   scale_fill_manual(values = c("Proto+Omi+" = "#648FFF",
#                                "Proto+Omi-" = "#FFB000",
#                                "Proto-Omi+" = "#DC267F"))+
#   theme_void()+
#   guides(fill = guide_legend(override.aes = list(size=3)))+
#   theme(legend.key.size = unit(0.1, 'cm'),
#         legend.text = element_text(size = 7),
#         legend.title = element_blank())
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_UMAP_ByProbeSpecificity.png"),width = 3.2, height = 3, units = "in", device = "png", dpi = 600)
# dev.off()
#####


#####
#Do a tabulated set on a per-donor basis, looking at proportional ratio of AM1/2 P+O+ : AM1/2 P+O at day 15
df$ClusterLabelAdj <- ifelse(df$ClusterLabel %in% c("AM1 (Activated)", "AM2 (Intermediate)"), "Activated", df$ClusterLabel)

stats <- df %>% filter(adj.ProtoOmi != "Proto-Omi+" & Timepoint == "Day 15") %>%
  group_by(OfficialBooster, Subject, adj.ProtoOmi, Timepoint, ClusterLabelAdj) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  filter(ClusterLabelAdj == "Activated") %>%
  ungroup()%>%
  group_by(OfficialBooster, Subject, Timepoint, ClusterLabelAdj) %>%
  arrange(adj.ProtoOmi) %>%
  mutate(PropRatio = Proportion[1] / Proportion[2],
         Booster = case_when(OfficialBooster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                             OfficialBooster == "Prototype + Omicron BA.1 mRNA" ~ "Prototype + BA.1",
                             OfficialBooster == "Prototype mRNA" ~ "Prototype"))

ggplot(stats[stats$adj.ProtoOmi != "Proto+Omi-",], aes(x = Booster, y = PropRatio))+
  geom_point(aes(fill = Booster), shape =21, size = 1.2, stroke = 0.8, position = position_jitter(width = 0.1))+
  ylab("Cross-Reactive / Prototype-Specific")+
  ggtitle("Activation by Probe Specificity")+
  scale_fill_manual(values = c("Omicron BA.1" = "#7C1D6f", 
                               "Prototype + BA.1" = "#DC3977",
                               "Prototype" = "#045275"))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        title = element_text(size = 8),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_ProportionalRatiosOfActivated_CrossVsProto.png"),width = 2.2, height = 2.5, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Heatmap 
#Part 1:
#Heatmap showing IGC usage 
stats <- naiveDF %>% filter(!is.na(c_call) & Timepoint == "Day 180") %>%
  group_by(ClusterLabel, c_call) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n)) %>%
  ungroup() %>%
  complete(ClusterLabel, c_call, fill = list(Proportion = 0))

ggplot(stats, aes(y=c_call, x=ClusterLabel, alpha=Proportion))+
  geom_tile(color="white", linewidth=1, fill="dodgerblue4")+
  #geom_tile(color="white", linewidth=1, aes(fill=ClusterLabel))+
  scale_alpha(range = c(0.1, 1.4), limits = c(0, 0.8))+
  #scale_fill_manual(values = phenoColors)+
  coord_equal()+
  #scale_y_discrete(limits=rev(c("Week 0","Week 1","Week 4", "Week 16", "Week 17", "Week 18", "Week 20", "Week 28", "Week 40")))+
  scale_y_discrete(limits = c("IGHM", "IGHD", "IGHG1", "IGHG3", "IGHG2", "IGHG4", "IGHA1", "IGHA2"))+
  scale_x_discrete(limits = rev(c("AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Atypical IgA Memory", "Resting IgA Memory", "Resting IgG Memory", "Unclear", "Naive")), position = "top")+
  #labs(y="Cluster")+
  #scale_fill_brewer(palette="Blues")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0, hjust=0, size=6), 
        axis.text.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
        legend.position = "left")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_FCProportionalUsage.png"),width = 3.4, height = 2.5, units = "in", device = "png", dpi = 600)
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
                filter(marker %in% c("IgA", "IgG", "IgM", "CD21", "CD27", "CD11c", "CD62L", "CD85j", "FCRL5", "CD72", "CD95", "CD45RB")) %>%
                mutate(marker = factor(marker, levels = c("IgG", "IgA", "IgM", "CD21", "CD27", "CD62L", "CD85j", "CD11c", "FCRL5", "CD72", "CD95", "CD45RB")))

ggplot(dfExtended, aes(x= ClusterLabel, y = expression))+
  geom_violin(aes(fill = ClusterLabel), linewidth = 0.2)+
  scale_fill_manual(values = phenoColors)+
  scale_y_continuous(position = "right", n.breaks = 3)+
  scale_x_discrete(limits = rev(c("AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)", "Atypical IgA Memory", "Resting IgA Memory", "Resting IgG Memory", "Unclear", "Naive")))+
  facet_grid(rows = vars(marker), switch="y", scales= "free_y")+
  annotate("rect", xmin=-Inf, xmax = Inf, ymin= -Inf, ymax = Inf, fill = NA, color = "black")+
  geom_vline(xintercept = -Inf, linewidth = 0.8)+
  theme_classic()+
  theme(legend.text = element_text(size = 7),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing.y = unit(0, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_ProteinViolinPlot.png"),width = 2.55, height = 3.3, units = "in", device = "png", dpi = 600)

###making gene expression plot
# allmarkers <- FindAllMarkers(seuObj, assay = "RNA", min.pct = .05, logfc.threshold = 0.3)
# subsetMarkers <- allmarkers %>% filter(abs(avg_log2FC) > 4 & allmarkers$p_val > 0 & allmarkers$p_val < 0.05)
# DotPlot(seuObj, features = rownames(subsetMarkers), split.by="ClusterLabel", assay = "RNA", cols = c("#D53E4F", "#F46D43","#FDAE61", "#FEE08B","#E6F598",
#                                                                                                 "#ABDDA4","#66C2A5","#3288BD"))

#####


######
###New plots
######
#Look at distribution of phenotypes by probe population
stats <- df %>% group_by(adj.ProtoOmi, Timepoint, ClusterLabel) %>%
          mutate(ClusterLabel = factor(ClusterLabel, levels = c("AM1 (Activated)",
                                                                 "AM2 (Intermediate)",
                                                                 "AM3 (Atypical)",
                                                                 "Atypical IgA Memory",
                                                                 "Resting IgA Memory",
                                                                 "Resting IgG Memory",
                                                                 "Unclear"))) %>%
          filter(Infection == "N") %>%
          summarize( n = n()) %>%
          mutate(Proportion = n / sum(n))

ggplot(stats, aes(x = Timepoint, y = Proportion, alluvium = ClusterLabel, fill = ClusterLabel, label = ClusterLabel, stratum = ClusterLabel))+
  geom_flow()+
  geom_stratum(linewidth = 0.4)+
  scale_fill_manual(values = phenoColors)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  facet_grid(cols = vars(adj.ProtoOmi), labeller = label_wrap_gen(width = 15))+
  theme(legend.text = element_text(size=6),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "lines"),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=0),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 7),
        #legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size=7))+
  guides(shape = guide_legend(override.aes = list(size=0.4)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_AlluvialPlot_ByProbeSpec.png"),width = 4.5, height = 2.1, units = "in", device = "png", dpi = 600)

#umap over time
naiveDF$Timepoint <- factor(naiveDF$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = ClusterLabel))+
  geom_point(size = 0.8, shape=21, stroke = 0.1)+
  scale_fill_manual(values = phenoColors)+
  facet_grid(cols = vars(Timepoint))+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size=3)))+
  theme(legend.key.size = unit(0.1, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_blank())
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4a_UMAP_overtime.png"),width = 4.3, height = 1.3, units = "in", device = "png", dpi = 600)
dev.off()

#probe specificity by population
ggplot(df, aes(x = ClusterLabel, fill= adj.ProtoOmi))+
  geom_bar(position = "fill")+
  scale_x_discrete(limits = c("AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)",
                              "Resting IgG Memory", "Unclear",
                              "Atypical IgA Memory", "Resting IgA Memory"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#8fc0a9",
                               "Proto-Omi+" = "#faf3dd",
                               "Proto+Omi-" = "#c8d5b9"))+
  ylab("Proportion")+
  xlab("Cluster")+
  ggtitle("Breadth by Cluster")+
  theme_classic()+
  theme(axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 8),
        axis.title.y = element_text(size = 8),
        title = element_text(size = 8))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProbeSpecificityByPopulation.png"),width = 4.9, height = 2.7, units = "in", device = "png", dpi = 800)


####sanity checking myself
#key
# stats <- df %>%
#         mutate(ClusterLabel = case_when(ClusterLabel == "Resting IgG Memory 2" ~ "Resting IgG Memory",
#                                         TRUE ~ ClusterLabel)) %>%
#         #filter(!(ClusterLabel %in% c("Resting IgA Memory", "Atypical IgA Memory") & c_call %in% c("IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"))) %>%
#         group_by(Booster, Subject, InfectionRange, Timepoint, ClusterLabel) %>%
#         summarize(n = n()) %>%
#         mutate(Proportion = n / sum(n))
# 
# ggplot(stats[stats$ClusterLabel == "Atypical IgA Memory",], aes(x = Timepoint, y = n))+
#   geom_point(shape = 21, aes(fill = Booster), stroke = 0.3)+
#   geom_line(aes(group = Subject))+
#   ggtitle("Activated IgA Over Time")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
#   theme_classic()+
#   theme(legend.key.size = unit(0.1, 'cm'),
#         legend.text = element_text(size = 6),
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ActivatedIgA_OverTime.png"),width = 6, height = 5.3, units = "in", device = "png", dpi = 600)
# 
# ggplot(stats[stats$ClusterLabel == "Resting IgA Memory",], aes(x = Timepoint, y = n))+
#   geom_point(shape = 21, aes(fill = Booster), stroke = 0.3)+
#   geom_line(aes(group = Subject))+
#   ggtitle("Resting IgA Over Time")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
#   theme_classic()+
#   theme(legend.key.size = unit(0.1, 'cm'),
#         legend.text = element_text(size = 6),
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "RestingIgA_OverTime.png"),width = 6, height = 5.3, units = "in", device = "png", dpi = 600)
# 
# ###
# stats <- df %>%
#   mutate(ClusterLabel = case_when(ClusterLabel == "Resting IgG Memory 2" ~ "Resting IgG Memory",
#                                   TRUE ~ ClusterLabel)) %>%
#   filter(!(ClusterLabel %in% c("Resting IgA Memory", "Atypical IgA Memory") & c_call %in% c("IGHD", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM"))) %>%
#   group_by(Booster, Subject, InfectionRange, Timepoint, ClusterLabel) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n))
# 
# ggplot(stats[stats$ClusterLabel == "Atypical IgA Memory",], aes(x = Timepoint, y = n))+
#   geom_point(shape = 21, aes(fill = Booster), stroke = 0.3)+
#   geom_line(aes(group = Subject))+
#   ggtitle("Activated IgA Over Time (True IgA)")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
#   theme_classic()+
#   theme(legend.key.size = unit(0.1, 'cm'),
#         legend.text = element_text(size = 6),
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ActivatedIgA_OverTime_cleaned.png"),width = 6, height = 5.3, units = "in", device = "png", dpi = 600)
# 
# ggplot(stats[stats$ClusterLabel == "Resting IgA Memory",], aes(x = Timepoint, y = n))+
#   geom_point(shape = 21, aes(fill = Booster), stroke = 0.3)+
#   geom_line(aes(group = Subject))+
#   ggtitle("Resting IgA Over Time (True IgA)")+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
#   facet_grid(cols = vars(Booster), rows = vars(InfectionRange))+
#   theme_classic()+
#   theme(legend.key.size = unit(0.1, 'cm'),
#         legend.text = element_text(size = 6),
#         legend.title = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "RestingIgA_OverTime_cleaned.png"),width = 6, height = 5.3, units = "in", device = "png", dpi = 600)
# 
# #show that they upregulate LGALS1, associated with weird memory subset (anergy?)
# VlnPlot(seuObj, features = "LGALS1", group.by = "ClusterLabel", pt.size = 0.01, assay="RNA")+
#   scale_x_discrete(limits = c("Atypical IgA Memory", "Resting IgA Memory",
#                               "AM1 (Activated)", "AM2 (Intermediate)", "AM3 (Atypical)",
#                               "Resting IgG Memory", "Resting IgG Memory 2", "Unclear"))+
#   theme(legend.position = "none")
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "LGALS1_expression.png"),width = 6, height = 4, units = "in", device = "png", dpi = 600)
