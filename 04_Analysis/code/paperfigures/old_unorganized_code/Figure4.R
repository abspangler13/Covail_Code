library(Seurat)
library(dplyr)
library(tidyseurat)
library(ggplot2)
library(ggalluvial)
library(here)
library(stringr)
library(writexl)

#load in the data
seuObjNaive <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds")) %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"),
         ClusterLabel = factor(ClusterLabel, levels = c("Atypical", "Acute Activated", "Intermediate", "Resting IgG", "Resting IgA", "Plasmablast-like", "Naive")),
         Booster = str_replace(Booster, "Omicron", "BA.1"),
         OfficialBooster = factor(OfficialBooster, levels = c("Prototype mRNA", "Prototype + Omicron BA.1 mRNA", "Omicron BA.1 mRNA")),
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")))

seuObj <- seuObjNaive %>% filter(ClusterLabel != "Naive", ClusterLabel != "Plasmablast-like")

df <- seuObj@meta.data %>% filter(adj.ProtoOmi != "Proto-Omi+")

df$OfficialInfection <- ifelse(df$Infection == "Y", "Infected", "Uninfected")
df$OfficialInfection <- factor(df$OfficialInfection, levels = c("Uninfected",
                                                                "Infected"))

df$InfectionRange <- ifelse(is.na(df$InfectionRange), "Uninfected", df$InfectionRange)

#set colors
allColors3 <- c("Prototype" = "#FBB042",
                     "BA.1 And Prototype" = "#1D75BC",
                     "BA.1" = "#2AB673")

allColors <- c("Omicron BA.1 mRNA" = "#2AB673", 
               "Prototype + Omicron BA.1 mRNA" = "#1D75BC",
               "Prototype mRNA" = "#FBB042")

allColors2 <- c("Omicron" = "#2AB673", 
                "Omicron And Prototype" = "#1D75BC",
                "Prototype" = "#FBB042")

shortColors <- c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 #"Acute Activated" = "#F46D43",
                 "Acute Activated" = "#f08665",
                 "Intermediate" = "#E6F598",
                 "Resting IgG" = "limegreen",
                 "Resting IgA" = "#3288BD",
                 "Plasmablast-like" = "#6f2da8",
                 "Naive" = "white")

naiveDF <- seuObjNaive@meta.data
embeddings <- as.data.frame(Embeddings(seuObjNaive,reduction = "harmony.wnn.umap"))
naiveDF$UMAP_1 <- embeddings$harmonywnnUMAP_1[match(rownames(embeddings), rownames(naiveDF))]
naiveDF$UMAP_2 <- embeddings$harmonywnnUMAP_2[match(rownames(embeddings), rownames(naiveDF))]

#create table
table(naiveDF$Booster, naiveDF$Timepoint)
######

#####
#Figure 4: showing clustering results of all cells
#let's plot this out- first
ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = ClusterLabel))+
  geom_point(size = 0.9, shape=21, stroke = 0.1)+
  scale_fill_manual(values = shortColors)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  theme(legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(0.1, "line"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4a_UMAP.png"),width = 2.9, height = 2.7, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4a_UMAP.svg"),width = 2.9, height = 2.7, units = "in")
dev.off()
#####

#####
#Figure 3b: showing phenotypic changes over time but split by specificity
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4b_AlluvialPlotOfPopulations_splitbyspecificity.png"),width = 4.2, height = 2.6, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4b_AlluvialPlotOfPopulations_splitbyspecificity.svg"),width = 4.2, height = 2.6, units = "in", device = "svg")
dev.off()
#####

#####
#Make a general summary alluvial plot that includes all cells, all boosters, all phenotypes
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "GeneralAlluvialPlot.png"),width = 3.2, height = 2.2, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "GeneralAlluvialPlot.svg"),width = 3.2, height = 2.2, units = "in", device = "svg")
dev.off()
#####

#####
#Heatmap 
#Part 1:
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_FCProportionalUsage.png"),width = 2.7, height = 3.3, units = "in", device = "png", dpi = 1200)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_FCProportionalUsage.svg"),width = 2.7, height = 3.3, units = "in")
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_ProteinViolinPlot_reduced.png"),width = 4.2, height = 2.2, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "Figure4_ProteinViolinPlot_reduced.svg"),width = 4.2, height = 2.2, units = "in")
dev.off()
#####

#####
#What proportion of activated cells post-boost are omicron-binding
# stats <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
#   mutate(RBDBinding = ifelse(adj.ProtoOmi == "Proto+Omi-", "Prototype-Specific", "Omicron-Binding"),
#          Booster = str_replace(Booster, "Omicron", "BA.1")) %>%
#   filter(Timepoint == "Day 15") %>%
#   group_by(Booster, Subject, ClusterLabel,RBDBinding) %>%
#   summarize(n = n()) %>%
#   mutate(Proportion = n / sum(n)) %>%
#   filter(RBDBinding == "Omicron-Binding")
# 
# stats %>% filter(ClusterLabel == "Acute Activated") %>%
#   ggplot(aes(x = Booster, y = Proportion))+
#   geom_point(shape = 21, aes(fill = Booster), position = position_jitter(width=0.1))+
#   scale_fill_manual(values = c("BA.1" = "#7C1D6f", 
#                                "BA.1 And Prototype" = "#DC3977",
#                                "Prototype" = "#076894"))+
#   scale_x_discrete(limits = c("Prototype", "BA.1 And Prototype", "BA.1"))+
#   ylim(0,1.05)+
#   ylab("Proportion Omicron+ at Day 15")+
#   ggtitle("Acute Activated")+
#   theme_classic()+
#   theme(legend.position = "none",
#         axis.text.x = element_text(hjust =1, vjust=1, angle = 45, size = 7),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 7),
#         axis.title.y = element_text(size = 7),
#         title = element_text(size = 7))
# ggsave(filename = here::here("04_Analysis", "plots","paperfigures", "Figure 4", "AcuteActivated_OmicronPositive.png"), width = 1.3, height = 1.9)
# ggsave(filename = here::here("04_Analysis", "plots","paperfigures", "Figure 4", "AcuteActivated_OmicronPositive.svg"), width = 1.3, height = 1.9)
# 
# #intermediate
# stats %>% filter(ClusterLabel == "Intermediate") %>%
#   ggplot(aes(x = Booster, y = Proportion))+
#   geom_point(shape = 21, aes(fill = Booster), position = position_jitter(width=0.1))+
#   scale_fill_manual(values = c("BA.1" = "#7C1D6f", 
#                                "BA.1 And Prototype" = "#DC3977",
#                                "Prototype" = "#076894"))+
#   ylim(0,1.05)+
#   ylab("Proportion Omicron+ at Day 15")+
#   scale_x_discrete(limits = c("Prototype", "BA.1 And Prototype", "BA.1"))+
#   ggtitle("Intermediate")+
#   theme_classic()+
#   theme(legend.position = "none",
#         axis.text.x = element_text(hjust =1, vjust=1, angle = 45, size = 7),
#         axis.title.x = element_blank(),
#         axis.text.y = element_text(size = 7),
#         axis.title.y = element_text(size = 7),
#         title = element_text(size = 7))
# ggsave(filename = here::here("04_Analysis", "plots","paperfigures", "Figure 4", "Intermediate_OmicronPositive.png"), width = 1.3, height = 1.9)
# ggsave(filename = here::here("04_Analysis", "plots","paperfigures", "Figure 4", "Intermediate_OmicronPositive.svg"), width = 1.3, height = 1.9)
# 
# #write a csv
# stats %>% filter(ClusterLabel %in% c("Acute Activated", "Intermediate")) %>%
#   write.csv(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 4", "proportion_activated_omicronbinding.csv"))
# ###
# 
# ####
# #look at fold changes in cross-reactivity per population
# stats <- naiveDF %>% filter(Timepoint %in% c("Day 0", "Day 15"), !(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) %>%
#   mutate(RBDBinding = ifelse(adj.ProtoOmi == "Proto+Omi-", "Prototype-Specific", "Omicron-Binding")) %>%
#   group_by(Booster, Subject, ClusterLabel, Timepoint, RBDBinding) %>%
#   summarize(n = n()) %>%
#   mutate(Prop = n / sum(n)) %>%
#   filter(RBDBinding == "Omicron-Binding") %>%
#   group_by(Booster, Subject, ClusterLabel) %>%
#   arrange(Timepoint) %>%
#   mutate(Fold = Prop / Prop[1]) %>%
#   filter(Timepoint == "Day 15", ClusterLabel %in% c("Acute Activated", "Intermediate"))
# 
# #acute activated
# stats %>% filter(ClusterLabel == "Acute Activated") %>%
#   ggplot(aes(x = Booster, y = Fold))+
#   geom_point(aes(fill = Booster, group = Booster), shape = 21, position = position_jitter(width = 0.15), size = 1.3)+
#   geom_hline(yintercept = 1, linetype = 2)+
#   ylim(0,3.5)+
#   ylab("Fold Change")+
#   ggtitle("Cross-Reactive")+
#   scale_fill_manual(values = allColors3)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         axis.text.x = element_text(angle = 45, hjust =1 , vjust =1),
#         legend.position = "none",
#         axis.title.x = element_blank())
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldBA1Increases_AcuteActivated.png"), width= 1.2, height =1.9)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldBA1Increases_AcuteActivated.svg"), width= 1.2, height =1.9)
# 
# #intermediate
# stats %>% filter(ClusterLabel == "Intermediate") %>%
#   ggplot(aes(x = Booster, y = Fold))+
#   geom_point(aes(fill = Booster, group = Booster), shape = 21, position = position_jitter(width = 0.15), size = 1.3)+
#   geom_hline(yintercept = 1, linetype = 2)+
#   ylim(0,3.5)+
#   ylab("Fold Change")+
#   ggtitle("Cross-Reactive")+
#   scale_fill_manual(values = allColors3)+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         axis.text.x = element_text(angle = 45, hjust =1 , vjust =1),
#         legend.position = "none",
#         axis.title.x = element_blank())
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldBA1Increases_Intermediate.png"), width= 1.2, height =1.9)
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldBA1Increases_Intermediate.svg"), width= 1.2, height =1.9)
# 
# #write stats sheet
# stats %>% select(!c(Timepoint, RBDBinding, n, Prop)) %>%
#   write.csv(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 4", "FC_OmicronBinding_PerCluster.csv"))
# 
# dev.off()
# #####

#####
#look at how cross-reactive or proto-spec. activated responses proceed per donor
#Combine intermediate and activated?
calcs <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), adj.ProtoOmi != "Proto-Omi+") %>%
  mutate(ClusterLabel = case_when(ClusterLabel %in% c("Intermediate", "Acute Activated") ~ "Activated",
                                  TRUE ~ ClusterLabel)) %>%
  mutate(ClusterSpec = paste0(ClusterLabel,"_", adj.ProtoOmi)) %>%
  group_by(Booster, Subject, Timepoint, ClusterSpec) %>%
  summarize(n = n()) %>%
  group_by(Booster)%>%
  complete(Subject, Timepoint, ClusterSpec, fill = list(n = 0)) %>%
  group_by(Booster, Subject, Timepoint)%>%
  mutate(Prop = n / sum(n)) %>%
  filter(!is.na(ClusterSpec)) %>%
  group_by(Booster, Timepoint, ClusterSpec) %>%
  summarize(n = n(),
            mean = mean(Prop),
            sd = sd(Prop)) %>%
  mutate(ClusterLabel = str_extract(ClusterSpec, "Activated"),
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Cross-Reactive", "Prototype-Specific"),
         se = sd / sqrt(n)) %>%
  filter(!is.na(ClusterLabel), Timepoint %in% c("Day 0", "Day 15"))

calcs %>% filter(ClusterLabel == "Activated") %>%
  ggplot(aes(x = Timepoint, y = mean))+
  geom_line(aes(group = adj.ProtoOmi, color = adj.ProtoOmi), alpha = 0.9)+
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se, color = adj.ProtoOmi), width = 0.3)+
  geom_point(shape =21, aes(fill = adj.ProtoOmi), size = 1)+
  scale_x_discrete(limits= c("Day 0", "Day 15"))+
  ylab("Proportion of Total Cells")+
  ggtitle("Activated")+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Prototype-Specific" = "lightgray",
                                "Cross-Reactive" = "#343148ff"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 7),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        plot.title = element_text(hjust = 0.5, size = 6))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "AcuteActivated_CrossvsProto_combinedclusters.png"), width = 3.1, height = 1.9, dpi = 800)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "AcuteActivated_CrossvsProto_combinedclusters.svg"), width = 3.1, height = 1.9)

#write xlsx
calcs <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), adj.ProtoOmi != "Proto-Omi+") %>%
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
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Cross-Reactive", "Prototype-Specific")) %>%
  filter(!is.na(ClusterLabel)) %>%
  select(Booster, Subject, adj.ProtoOmi, Timepoint, Prop) %>%
  pivot_wider(names_from = Timepoint, values_from = Prop) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure 4", "ClustersCombined_ProportionActivatedCells.xlsx"))


######now calculate FC/Delta
calcs <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), adj.ProtoOmi != "Proto-Omi+") %>%
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
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Cross-Reactive", "Prototype-Specific")) %>%
  filter(!is.na(ClusterLabel)) %>%
  group_by(Booster, Subject, ClusterSpec) %>%
  arrange(Timepoint) %>%
  mutate(Delta = Prop - Prop[1], FC = Prop / Prop[1]) %>% filter(Timepoint == "Day 15", !is.infinite(FC))

#plot
calcs %>% filter(ClusterLabel == "Activated") %>%
  ggplot(aes(x = Booster, y = FC))+
  geom_point(shape = 21, aes(fill = Booster), size = 1, position=  position_jitter(width = 0.2))+
  scale_fill_manual(values = c(allColors3))+
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3)+
  facet_grid(cols = vars(adj.ProtoOmi))+
  ylab("Fold Change")+
  ylim(0, 10)+
  theme_classic()+
  theme(text = element_text(size = 7),
        legend.position = "none",
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, size = 5, hjust = 1, vjust = 1),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldChange_ProportionActivated_combinedclusters.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "FoldChange_ProportionActivated_combinedclusters.svg"), width = 1.6, height = 1.8)


#write stats sheet
calcs %>% filter(ClusterLabel == "Activated")%>% ungroup() %>% select(Booster, Subject, adj.ProtoOmi, Delta, FC) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures","Figure 4", "AcuteActivated_Delta_bySpecificity_combinedclusters.xlsx"))
########

########
#Alternative plot to above: show differences per-donor as proportion of total response
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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_PrototypeBoost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_PrototypeBoost.svg"), width = 1.6, height = 1.8)

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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_BivalentBoost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_BivalentBoost.svg"), width = 1.6, height = 1.8)

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
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_BA1Boost.png"), width = 1.6, height = 1.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "ProportionTotalCells_Activated_BA1Boost.svg"), width = 1.6, height = 1.8)
######

######
#Alternative: plotting proportion of activated cells that are cross-reactive vs prototype-specific
calcs <- df %>% filter(Timepoint %in% c("Day 15"), ClusterLabel %in% c("Acute Activated", "Intermediate")) %>%
        group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
        summarize(n = n()) %>%
        mutate(Prop = n / sum(n),
               adj.ProtoOmi = ifelse(adj.ProtoOmi == "Proto+Omi+", "Cross-Reactive", "Prototype-Specific"))

#plot out the differences per group
#prototype
calcs %>%
  ggplot(aes(x = adj.ProtoOmi, y = Prop))+
  geom_line(aes(group = Subject))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21)+
  ylab("Proportion of Activated Cells, Day 15")+
  ggtitle("Prototype Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(Booster))+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Prototype-Specific" = "lightgray",
                                "Cross-Reactive" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        plot.title = element_text(hjust = 0.5))

###try whole number
calcs %>%
  ggplot(aes(x = adj.ProtoOmi, y = n))+
  geom_line(aes(group = Subject))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21)+
  ylab("Proportion of Activated Cells, Day 15")+
  ggtitle("Prototype Boost")+
  facet_grid(cols = vars(Booster))+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Prototype-Specific" = "lightgray",
                                "Cross-Reactive" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1),
        plot.title = element_text(hjust = 0.5))

#proto/BA1
calcs %>% filter(Booster == "BA.1 And Prototype") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21)+
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
        plot.title = element_text(hjust = 0.5))

#BA.1 vaccine
calcs %>% filter(Booster == "BA.1") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi))+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21)+
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
        plot.title = element_text(hjust = 0.5))
#######

#######
#Plot the relative phenotypes of cells at day 180 that are BA.1-specific
stats <- naiveDF %>% filter(Timepoint == "Day 180", Infection == "N", adj.ProtoOmi == "Proto-Omi+") %>% group_by(Booster, adj.ProtoOmi, ClusterLabel) %>%
        summarize(n = n()) %>% mutate(Prop = n / sum(n), sum = sum(n))

ggplot(stats, aes(x = Booster, y = Prop, fill = ClusterLabel))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  geom_text(mapping = aes(label = sum, x = Booster, y = 1.05), size = 3)+
  ylim(0,1.05)+
  scale_fill_manual(values = shortColors)+
  ylab("Proportion")+
  ggtitle("Prototype-BA.1+, Day 180")+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "BA1Positive_d180_Phenotypes.png"), width = 2.4, height = 2.4)

#plot phenotypes of BA.1+ b cells from donors with a particularly large BA.1+ increase by d180
check <- naiveDF %>% filter(Infection == "N", Timepoint == "Day 180") %>% group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>% summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>%
            filter(adj.ProtoOmi == "Proto-Omi+")

denovoDonors <- check$Subject[check$Proportion >= 0.1]

#do the same as above except for these donors only
stats <- naiveDF %>% filter(Subject %in% denovoDonors, adj.ProtoOmi == "Proto-Omi+", Timepoint == "Day 180") %>%
          group_by(Booster, Subject, ClusterLabel) %>% summarize(n = n()) %>% mutate(Prop = n / sum(n), sum = sum(n))

ggplot(stats, aes(x = Subject, y = Prop))+
  geom_bar(stat = "identity", position = "stack", aes(fill = ClusterLabel), color = "black")+
  geom_text(mapping = aes(label = sum, x = Subject, y = 1.05), size = 3)+
  ylim(0,1.05)+
  scale_fill_manual(values = shortColors)+
  facet_grid(cols = vars(Booster), scale = "free_x")+
  ylab("Proportion")+
  ggtitle("Prototype-BA.1+, Day 180")+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        legend.key.size = unit(0.8, "lines"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure 4", "BA1Positive_d180_Phenotypes_expanded_donors.png"), width = 4, height = 2.4)
