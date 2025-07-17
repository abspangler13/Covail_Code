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
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")),
         Boost = case_when(Booster == "BA.1 And Prototype" ~ "Prototype/BA.1",
                           TRUE ~ Booster),
         Boost = factor(Boost, levels = c("Prototype", "Prototype/BA.1", "BA.1")))

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

allColors4 <- c("BA.1" = "#2AB673", 
                "Prototype/BA.1" = "#1D75BC",
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

######
#plot citeseq specificities over time
stats <- naiveDF %>%
  group_by(Boost, Infection, InfectionRange, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Subject = as.character(Subject)) %>%
  filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180"))) #let's include infected donors but only at days 0 and 15

ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi+",], aes(x=Timepoint, y=Proportion, fill=Boost))+
  geom_line(alpha = 0.5, aes(group = Subject, color = Boost))+
  geom_point(shape=21, aes(fill=Boost), stroke = 0.5)+
  ylab("Proportion of All Cells")+
  ggtitle("Prototype+BA.1+")+
  facet_grid(cols=vars(Boost), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors4)+
  scale_color_manual(values = allColors4)+
  theme_classic() +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "CITESeq_CrossReactives.svg"),width = 5.2, height = 3.5)
dev.off()

#prototype specific
ggplot(stats[stats$adj.ProtoOmi == "Proto+Omi-",], aes(x=Timepoint, y=Proportion, fill=Boost))+
  geom_line(alpha = 0.5, aes(group = Subject, color = Boost))+
  geom_point(shape=21, aes(fill=Boost), stroke = 0.5)+
  ylab("Proportion of All Cells")+
  ggtitle("Prototype+BA.1-")+
  facet_grid(cols=vars(Boost), labeller = label_wrap_gen(15))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors4)+
  scale_color_manual(values = allColors4)+
  theme_classic() +
  theme(text = element_text(size = 16),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "CITESeq_PrototypeSpecifics.svg"),width = 5.2, height = 3.5)
dev.off()
######

######
#plot UMAP
ggplot(naiveDF, aes(x = UMAP_1, y= UMAP_2, fill = ClusterLabel))+
  geom_point(size = 0.9, shape=21, stroke = 0.015)+
  scale_fill_manual(values = shortColors)+
  theme_void()+
  guides(fill = guide_legend(override.aes = list(size = 4)))+
  theme(legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(0.1, "line"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "UMAP.svg"),width = 4.5, height = 3.8, units = "in")
dev.off()
######

######
#Alluvial plot
stats <- naiveDF %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), adj.ProtoOmi != "Proto-Omi+") %>%
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
  theme(legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = "right",
        legend.key.size = unit(1, "lines"),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=12,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size = 12),
        strip.background = element_blank())+
  guides(shape = guide_legend(override.aes = list(size=0.9)))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "PhenotypicAlluvium.svg"),width = 5, height = 3.5, units = "in")
dev.off()

######

######
#Show activation over time
#Combine intermediate and activated?
calcs <- df %>% filter(!(Infection == "Y" & Timepoint %in% c("Day 90", "Day 180")), adj.ProtoOmi != "Proto-Omi+") %>%
  mutate(ClusterLabel = case_when(ClusterLabel %in% c("Intermediate", "Acute Activated") ~ "Activated",
                                  TRUE ~ ClusterLabel)) %>%
  mutate(ClusterSpec = paste0(ClusterLabel,"_", adj.ProtoOmi)) %>%
  group_by(Boost, Subject, Timepoint, ClusterSpec) %>%
  summarize(n = n()) %>%
  group_by(Boost)%>%
  complete(Subject, Timepoint, ClusterSpec, fill = list(n = 0)) %>%
  group_by(Boost, Subject, Timepoint)%>%
  mutate(Prop = n / sum(n)) %>%
  filter(!is.na(ClusterSpec)) %>%
  group_by(Boost, Timepoint, ClusterSpec) %>%
  summarize(n = n(),
            mean = mean(Prop),
            sd = sd(Prop)) %>%
  mutate(ClusterLabel = str_extract(ClusterSpec, "Activated"),
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Cross-Reactive", "Prototype-Specific"),
         se = sd / sqrt(n)) %>%
  filter(!is.na(ClusterLabel), Timepoint %in% c("Day 0", "Day 15"))

calcs %>% filter(ClusterLabel == "Activated") %>%
  ggplot(aes(x = Timepoint, y = mean))+
  geom_line(aes(group = adj.ProtoOmi, color = adj.ProtoOmi), alpha = 0.9, linewidth = 0.8)+
  geom_errorbar(aes(ymax = mean + se, ymin = mean - se, color = adj.ProtoOmi), width = 0.3)+
  geom_point(shape =21, aes(fill = adj.ProtoOmi), size = 1.5)+
  scale_x_discrete(limits= c("Day 0", "Day 15"))+
  ylab("Proportion of Total Cells")+
  ggtitle("Activated Cells")+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  scale_color_manual(values = c("Prototype-Specific" = "lightgray",
                                "Cross-Reactive" = "#343148ff"))+
  facet_grid(cols = vars(Boost))+
  theme_classic()+
  theme(text = element_text(size = 14),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 10))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "Activated_OverTime_CombinedClusters.svg"), width = 6.7, height = 3.6)
######

######
#Alternative plot
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
         adj.ProtoOmi = ifelse(str_detect(ClusterSpec, "Proto\\+Omi\\+"), "Activated Prototype+BA.1+", "Activated Prototype+BA.1-"),
         adj.ProtoOmi = factor(adj.ProtoOmi, levels = c("Activated Prototype+BA.1+", "Activated Prototype+BA.1-"))) %>%
  filter(!is.na(ClusterLabel))

#plot out the differences per group
#prototype
calcs %>% filter(Booster == "Prototype") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi), linewidth = 0.8)+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 2)+
  ylab("Proportion of Total Cells")+
  ggtitle("Prototype Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(adj.ProtoOmi), labeller = label_wrap_gen(10))+
  scale_fill_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                               "Activated Prototype+BA.1+" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                                "Activated Prototype+BA.1+" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 15, color= "black"),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "ProportionTotalCells_Activated_perdonor_prototype.svg"), width = 3.6, height = 3.8)

#prototype/ba.1
calcs %>% filter(Booster == "BA.1 And Prototype") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi), linewidth = 0.8)+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 2)+
  ylab("Proportion of Total Cells")+
  ggtitle("Prototype/BA.1 Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(adj.ProtoOmi), labeller = label_wrap_gen(10))+
  scale_fill_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                               "Activated Prototype+BA.1+" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                                "Activated Prototype+BA.1+" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 15, color= "black"),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "ProportionTotalCells_Activated_perdonor_bivalent.svg"), width = 3.6, height = 3.8)

#ba.1
calcs %>% filter(Booster == "BA.1") %>%
  ggplot(aes(x = Timepoint, y = Prop))+
  geom_line(aes(group = Subject, color = adj.ProtoOmi), linewidth = 0.8)+
  geom_point(aes(fill = adj.ProtoOmi), shape = 21, size = 2)+
  ylab("Proportion of Total Cells")+
  ggtitle("BA.1 Boost")+
  ylim(0,1)+
  facet_grid(cols = vars(adj.ProtoOmi), labeller = label_wrap_gen(10))+
  scale_fill_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                               "Activated Prototype+BA.1+" = "#343148ff"))+
  scale_color_manual(values = c("Activated Prototype+BA.1-" = "lightgray",
                                "Activated Prototype+BA.1+" = "#343148ff"))+
  theme_classic()+
  theme(text = element_text(size = 15, color= "black"),
        strip.background = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust =1, color = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "ProportionTotalCells_Activated_perdonor_ba1.svg"), width = 3.6, height = 3.8)

######