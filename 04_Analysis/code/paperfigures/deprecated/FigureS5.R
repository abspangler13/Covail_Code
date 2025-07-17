library(Seurat)
library(dplyr)
library(tidyseurat)
library(ggplot2)
library(ggalluvial)
library(here)
library(stringr)
library(writexl)
library(readxl)
library(shazam)

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
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#076894")

allColors2 <- c("Omicron" = "#7C1D6f", 
                "Omicron And Prototype" = "#DC3977",
                "Prototype" = "#076894")

allColors3 <- c("BA.1" = "#7C1D6f", 
                "BA.1 And Prototype" = "#DC3977",
                "Prototype" = "#076894")

specColors <- c("Proto+Omi+" = "#386e72",
                "Proto+Omi-" = "#95C5C8",
                "Proto-Omi+" = "#F0C0AA")

shortColors <- c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                 #"Acute Activated" = "#F46D43",
                 "Acute Activated" = "#f08665",
                 "Intermediate" = "#E6F598",
                 "Resting IgG" = "limegreen",
                 "Resting IgA" = "#3288BD",
                 "Plasmablast-like" = "#6f2da8",
                 "Naive" = "white")

naiveDF <- seuObjNaive@meta.data

#read in evo data
evolving <- read.csv(here::here("01_raw-data", "Evo_dat_Timepoint_uniform.csv"))

#read in RATPIg-citeseq data
citeRat <- read_xlsx(here::here("01_raw-data", "CombinedMSD_Neut_Smartseq_OverlappingCITESeq.xlsx"))
######

######
#First, let's separate out significant lineages
evolvingLins <- evolving$clone_id[evolving$sig & !is.na(evolving$sig)]

naiveDF$sub_clone <- paste0("s", naiveDF$Subject, "_", naiveDF$clone_id, "_1")

evolved <- naiveDF %>% filter(sub_clone %in% evolvingLins)

#plot specificity over time
calcs <- evolved %>%
          group_by(Booster, adj.ProtoOmi, InfectionRange, clone_subject_id, Timepoint) %>%
          summarize(meanPrototype = mean(`Proto-RBD-PE`),
                    meanBA1 = mean(`BA1-RBD-PE`),
                    meanXBB = mean(`XBB-RBD-no-fluor`))

ggplot(calcs, aes(x = Timepoint, y = meanBA1))+
  geom_line(aes(group = clone_subject_id, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), axes = "all_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MeanBA1Signal_evolvingLineages.png"), height = 4, width = 4)

ggplot(calcs, aes(x = Timepoint, y = meanPrototype))+
  geom_line(aes(group = clone_subject_id, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), axes = "all_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MeanPrototypeSignal_evolvingLineages.png"), height = 4, width = 4)

ggplot(calcs, aes(x = Timepoint, y = meanXBB))+
  geom_line(aes(group = clone_subject_id, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), axes = "all_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MeanXBBSignal_evolvingLineages.png"), height = 4, width = 4)

#let's try calculating the proportion of labels that are prototype-specific per timepoint
calcs <- evolved %>%
        group_by(Booster, adj.ProtoOmi, clone_subject_id, Timepoint, ProtoOmi) %>%
        summarize(n = n()) %>% mutate(Proportion = n /sum(n))

calcs %>% filter(ProtoOmi == "Proto+Omi+")%>%
ggplot(aes(x = Timepoint, y = Proportion))+
  geom_line(aes(group = clone_subject_id, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ylab("Proportion Cross-Reactive Per Lineage")+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), axes = "all_x")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
######

######
#q2: overlaps b/w citeseq and ratpig for ba.1-specific B cells
deNovo <- citeRat %>% filter(str_detect(BindingPopulation, "(P-)|(BA.1-specific)"))

denovoCite <- naiveDF %>% filter(CELL %in% citeRat$CELL[citeRat$clone_subgroup_id %in% deNovo$clone_subgroup_id])

denovoEv <- evolving %>% filter(clone_id %in% denovoCite$sub_clone)

#plot relative group origin
calcs <- denovoCite %>% group_by(Booster) %>% summarize(n = n()) %>% mutate(P = n / sum(n))

ggplot(calcs, aes(x = 0, y = P))+
  geom_bar(stat = "identity", position = "stack", aes(fill = Booster))+
  scale_fill_manual(values = allColors3)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  ylab("Proportion of De Novo Overlap")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        text =element_text(size = 6),
        axis.ticks.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "ProportionOfDenovo.png"), width = 2.2, height = 3)

#probe labels without corrections
calcs <- denovoCite %>% group_by(clone_subject_id, ProtoOmi) %>% summarize(n = n()) %>% mutate(P = n / sum(n))

ggplot(calcs, aes(x = clone_subject_id, y = P))+
  geom_bar(stat = "identity", position = "stack", aes(fill = ProtoOmi))+
  #scale_fill_manual(values = allColors3)+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  ylab("Proportion of De Novo Overlap")+
  xlab("Clonal Group")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        text =element_text(size = 6))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "ProportionOfDenovo_specificity.png"), width = 4, height = 3.5)

#plot BA.1 signal over time for each group
calcs <- denovoCite %>% group_by(Booster, clone_subject_id, Timepoint) %>% summarize(mean = mean(`BA1-RBD-PE`))

ggplot(calcs, aes(x = Timepoint, y = mean))+
  geom_point(shape =21, aes(fill = Booster))+
  geom_line(aes(color = Booster, group = clone_subject_id))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ylab("Mean BA.1 Signal")+
  scale_color_manual(values = allColors3)+
  scale_fill_manual(values = allColors3)+
  theme_classic()+
  theme(text =element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust =1 , vjust = 1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MeanBA1Signal_RATPIgOverlap.png"), width = 3, height = 2.3)
########

########
#BA.1 boost over time
calcs <- naiveDF %>% filter(Infection == "N") %>%
          group_by(Booster, Timepoint, adj.ProtoOmi) %>%
          summarize(n = n()) %>% #mutate(Proportion = n / sum(n)) %>% select(!n) %>%
          #pivot_wider(names_from = Timepoint, values_from = Proportion)
          pivot_wider(names_from = Timepoint, values_from = n)

#ba1 booster
ba1 <- calcs %>% filter(Booster == "BA.1") %>% ungroup() %>% select(!c(Booster,`Day 0`, `Day 15`))
ba1 <- as.data.frame(ba1)
rownames(ba1) <- ba1$adj.ProtoOmi
ba1 <- ba1 %>% select(!adj.ProtoOmi)
ba1 <- as.matrix(ba1)
chisq.test(ba1)

#bivalent
biv <- calcs %>% filter(Booster == "BA.1 And Prototype") %>% ungroup() %>% select(!c(Booster,`Day 0`, `Day 15`))
biv <- as.data.frame(biv)
rownames(biv) <- biv$adj.ProtoOmi
biv <- biv %>% select(!adj.ProtoOmi)
biv <- as.matrix(biv)
chisq.test(biv)

#prototype
pro <- calcs %>% filter(Booster == "Prototype") %>% ungroup() %>% select(!c(Booster,`Day 0`, `Day 15`))
pro <- as.data.frame(pro)
rownames(pro) <- pro$adj.ProtoOmi
pro <- pro %>% select(!adj.ProtoOmi)
pro <- as.matrix(pro)
chisq.test(pro)

#####plot barplot for days 90 and 180
calcs <- naiveDF %>% filter(Infection == "N", Timepoint %in% c("Day 90", "Day 180")) %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>% select(!n)

#plot
ggplot(calcs, aes(x = Timepoint, y = Proportion))+
  geom_bar(position = "stack", stat = "identity", aes(fill = adj.ProtoOmi))+
  scale_x_discrete(limits=  c("Day 90", "Day 180"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1 , vjust =1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AntigenSpecificityBarplots.png"), width = 4, height = 3)


#######plot whole numbers over time
calcs <- naiveDF %>% filter(Infection == "N") %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n())

ggplot(calcs, aes(x = Timepoint, y = n))+
  geom_line(aes(group = Subject, color = Booster))+
  geom_point(aes(fill = Booster), shape = 21)+
  scale_x_discrete(limits= c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  ylab("Cells")+
  facet_grid(rows = vars(adj.ProtoOmi), cols = vars(Booster), axes = "all_x", scale= "free_y")+
  theme_classic()+
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45 , hjust =1 , vjust=1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "SpecificityCounts.png"), width = 4, height = 5)

#test whole numbers on day 90-day180 for BA.1+ b cells
calcs <- naiveDF %>% filter(Infection == "N", Timepoint %in% c("Day 90", "Day 180"), adj.ProtoOmi == "Proto-Omi+") %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  group_by(Booster)%>%
  complete(Subject, Timepoint, adj.ProtoOmi, fill = list(n = 0))

ggplot(calcs, aes(x = Timepoint, y = n))+
  geom_line(aes(group = Subject, color = Booster))+
  geom_point(aes(fill = Booster), shape = 21)+
  scale_x_discrete(limits= c("Day 90", "Day 180"))+
  scale_fill_manual(values = allColors3)+
  scale_color_manual(values = allColors3)+
  ggtitle("Prototype-Omicron+")+
  ylab("Cells")+
  facet_grid(cols = vars(Booster), axes = "all_x", scale= "free_y")+
  theme_classic()+
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45 , hjust =1 , vjust=1),
        strip.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "SpecificityCounts_ba1only_days90180.png"), width = 3, height = 2)

#test it
calcs <- naiveDF %>% filter(Infection == "N", Timepoint %in% c("Day 90", "Day 180"), adj.ProtoOmi == "Proto-Omi+") %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  group_by(Booster)%>%
  complete(Subject, Timepoint, adj.ProtoOmi, fill = list(n = 0)) %>%
  pivot_wider(names_from = Timepoint, values_from = n)

results <- by(calcs, calcs$Booster, function(x)
  wilcox.test(x$`Day 90`, x$`Day 180`, paired = TRUE))#[c(1:length(unique(calcs$Booster)))]
results <- type.convert(as.data.frame(do.call(rbind, results)), as.is=TRUE)
########

########
#Q3b: d1 vs d15?
calcs <- naiveDF %>% filter(Infection == "N") %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% #mutate(Proportion = n / sum(n)) %>% select(!n) %>%
  #pivot_wider(names_from = Timepoint, values_from = Proportion)
  pivot_wider(names_from = Timepoint, values_from = n)

#ba1 booster
ba1 <- calcs %>% filter(Booster == "BA.1") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
ba1 <- as.data.frame(ba1)
rownames(ba1) <- ba1$adj.ProtoOmi
ba1 <- ba1 %>% select(!adj.ProtoOmi)
ba1 <- as.matrix(ba1)
chisq.test(ba1)

#bivalent
biv <- calcs %>% filter(Booster == "BA.1 And Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
biv <- as.data.frame(biv)
rownames(biv) <- biv$adj.ProtoOmi
biv <- biv %>% select(!adj.ProtoOmi)
biv <- as.matrix(biv)
chisq.test(biv)

#prototype
pro <- calcs %>% filter(Booster == "Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
pro <- as.data.frame(pro)
rownames(pro) <- pro$adj.ProtoOmi
pro <- pro %>% select(!adj.ProtoOmi)
pro <- as.matrix(pro)
chisq.test(pro)

#####plot barplot for days 90 and 180
calcs <- naiveDF %>% filter(Infection == "N", Timepoint %in% c("Day 0", "Day 15")) %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>% select(!n)

#plot
ggplot(calcs, aes(x = Timepoint, y = Proportion))+
  geom_bar(position = "stack", stat = "identity", aes(fill = adj.ProtoOmi))+
  scale_x_discrete(limits=  c("Day 0", "Day 15"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1 , vjust =1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AntigenSpecificityBarplots_d0d15.png"), width = 4, height = 3)

#################with infected donors
calcs <- naiveDF %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% #mutate(Proportion = n / sum(n)) %>% select(!n) %>%
  #pivot_wider(names_from = Timepoint, values_from = Proportion)
  pivot_wider(names_from = Timepoint, values_from = n)

#ba1 booster
ba1 <- calcs %>% filter(Booster == "BA.1") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
ba1 <- as.data.frame(ba1)
rownames(ba1) <- ba1$adj.ProtoOmi
ba1 <- ba1 %>% select(!adj.ProtoOmi)
ba1 <- as.matrix(ba1)
chisq.test(ba1)

#bivalent
biv <- calcs %>% filter(Booster == "BA.1 And Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
biv <- as.data.frame(biv)
rownames(biv) <- biv$adj.ProtoOmi
biv <- biv %>% select(!adj.ProtoOmi)
biv <- as.matrix(biv)
chisq.test(biv)

#prototype
pro <- calcs %>% filter(Booster == "Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`))
pro <- as.data.frame(pro)
rownames(pro) <- pro$adj.ProtoOmi
pro <- pro %>% select(!adj.ProtoOmi)
pro <- as.matrix(pro)
chisq.test(pro)

#####plot barplot for days 90 and 180
calcs <- naiveDF %>% filter(Timepoint %in% c("Day 0", "Day 15")) %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>% select(!n)

#plot
ggplot(calcs, aes(x = Timepoint, y = Proportion))+
  geom_bar(position = "stack", stat = "identity", aes(fill = adj.ProtoOmi))+
  scale_x_discrete(limits=  c("Day 0", "Day 15"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1 , vjust =1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AntigenSpecificityBarplots_d0d15_withinfecteddonors.png"), width = 4, height = 3)


#########All 4 timepoints, uninfected donors
calcs <- naiveDF %>% filter(Infection == "N") %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>% select(!n)

#plot
ggplot(calcs, aes(x = Timepoint, y = Proportion))+
  geom_bar(position = "stack", stat = "identity", aes(fill = adj.ProtoOmi))+
  scale_x_discrete(limits=  c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Proto+Omi+" = "#386e72",
                               "Proto+Omi-" = "#95C5C8",
                               "Proto-Omi+" = "#F0C0AA"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust =1 , vjust =1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AntigenSpecificityBarplots_alltimepoints_withinfecteddonors.png"), width = 5, height = 3)

#do chi squared
calcs <- naiveDF %>% filter(Infection == "N") %>%
  group_by(Booster, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>% #mutate(Proportion = n / sum(n)) %>% select(!n) %>%
  #pivot_wider(names_from = Timepoint, values_from = Proportion)
  pivot_wider(names_from = Timepoint, values_from = n)

#ba1 booster
ba1 <- calcs %>% filter(Booster == "BA.1") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`, `Day 90`, `Day 180`))
ba1 <- as.data.frame(ba1)
rownames(ba1) <- ba1$adj.ProtoOmi
ba1 <- ba1 %>% select(!adj.ProtoOmi)
ba1 <- as.matrix(ba1)
chisq.test(ba1)

#bivalent
biv <- calcs %>% filter(Booster == "BA.1 And Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`, `Day 90`, `Day 180`))
biv <- as.data.frame(biv)
rownames(biv) <- biv$adj.ProtoOmi
biv <- biv %>% select(!adj.ProtoOmi)
biv <- as.matrix(biv)
chisq.test(biv)

#prototype
pro <- calcs %>% filter(Booster == "Prototype") %>% ungroup() %>% select(c(adj.ProtoOmi,`Day 0`, `Day 15`,`Day 90`, `Day 180`))
pro <- as.data.frame(pro)
rownames(pro) <- pro$adj.ProtoOmi
pro <- pro %>% select(!adj.ProtoOmi)
pro <- as.matrix(pro)
chisq.test(pro)

########

########
#Q4: Do we see differences in singlet vs non-singlet populations at each timepoint?
nonSinglets <- unique(naiveDF$clone_subject_id[duplicated(naiveDF$clone_subject_id) | duplicated(naiveDF$clone_subject_id, fromLast=T)])
naiveDF$CloneStatus <- ifelse(naiveDF$clone_subject_id %in% nonSinglets, "Expanded", "Singlet")

#plot
naiveDF %>% filter(Infection == "N") %>%
ggplot(aes(x = Timepoint, y = mu_freq * 100, fill = CloneStatus))+
  geom_violin()+
  geom_boxplot(width = 0.2, outlier.size = 0, position = position_dodge(width= 0.9))+
  scale_fill_manual(values = c("Singlet" = "white", "Expanded"= "skyblue"))+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  ylab("% VH Mutation")+
  ylim(0,12.5)+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 6),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1 ,vjust = 1))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "VHMut_ExpandedVsSinglet.png"), width = 6, height = 2.3, dpi = 1200)

#write a file
naiveDF %>% filter(Infection == "N") %>% select(Booster, Timepoint, CloneStatus, mu_freq) %>% mutate(mu_freq = mu_freq*100) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S5", "SHM_ByCloneStatus.xlsx"))

####look at it per donor over time
calcs <- naiveDF %>% filter(Infection == "N") %>% group_by(Booster, Subject, Timepoint, CloneStatus) %>%
          summarize(mean = mean(mu_freq * 100)) %>% group_by(Booster,Timepoint, CloneStatus) %>%
          summarize(mean2 = mean(mean), n = n(), sd = sd(mean)) %>%
          mutate(se = sd / sqrt(n))

ggplot(calcs, aes(x = Timepoint, y =mean2))+
  geom_line(aes(group = CloneStatus, color = CloneStatus))+
  geom_errorbar(aes(ymax = mean2+se, ymin = mean2-se, color = CloneStatus))+
  geom_point(aes(fill = CloneStatus), shape =21)+
  ylab("Mean % VH Mutation")+
  scale_fill_manual(values = c("Singlet" = "gray90", "Expanded"= "skyblue"))+
  scale_color_manual(values = c("Singlet" = "gray90", "Expanded"= "skyblue"))+
  facet_grid(cols = vars(Booster))+
  scale_x_discrete(limits= c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MeanVHMut_ByCloneStatus_OverTime.png"), width = 4, height =2.5, dpi = 1200)
########

########
#BA.1-specific changes in SHM over time?
#per donor
calcs <- naiveDF %>% filter(adj.ProtoOmi == "Proto-Omi+", Infection == "N") %>%
          group_by(Booster, Subject, Timepoint) %>%
          mutate(mu_freq = mu_freq * 100) %>%
          summarize(median = median(mu_freq))

ggplot(calcs, aes(x= Timepoint, y = median, color = Booster, fill = Booster))+
  geom_line(aes(group = Subject), linewidth = 0.6, alpha = 0.7)+
  geom_point(shape = 21)+
  ylab("Median % VH Mutation")+
  scale_color_manual(values = allColors3)+
  scale_fill_manual(values = allColors3)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1 , vjust =1),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "MedianVHMut_BA1Specific_PerDonor.png"), width = 4, height =2.5, dpi = 1200)

#overall
naiveDF %>% filter(Infection == "N") %>%
ggplot(aes(x= Timepoint, y = mu_freq*100, fill = adj.ProtoOmi))+
  geom_violin()+
  geom_boxplot(width = 0.4)+
  ylim(0, 13)+
  scale_fill_manual(values= specColors)+
  ylab("% VH Mutation")+
  facet_grid(rows = vars(Booster), cols = vars(adj.ProtoOmi), axes = "all_x")+
  scale_x_discrete(limits=c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.key.size = unit(0.5, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "VHMut_BySpecificity.png"), width = 6, height =5, dpi = 1200)

########Lineplot alternative?
calcs <- naiveDF %>% filter(Infection == "N") %>%
            group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
            mutate(mu_freq = mu_freq * 100) %>%
            summarize(mean = mean(mu_freq)) %>%
            group_by(Booster, Timepoint, adj.ProtoOmi) %>%
            summarize(n = n(), mean2 = mean(mean), sd = sd(mean)) %>%
            mutate(se = sd / sqrt(n))
            
ggplot(calcs, aes(x = Timepoint, y = mean2))+
  geom_errorbar(aes(ymin = mean2-se, ymax=mean2+se, color = adj.ProtoOmi))+
  geom_line(aes(color = adj.ProtoOmi, group = adj.ProtoOmi), linewidth = 0.8, alpha = 0.7)+
  geom_point(shape = 21, aes(fill = adj.ProtoOmi))+
  facet_grid(cols = vars(Booster))+
  ylab("Mean % VH Mutation")+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_color_manual(values = specColors)+
  scale_fill_manual(values = specColors)+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
        legend.key.size = unit(0.9, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "VHMut_BySpecificity_meanperdonor.png"), width = 4, height =2.2, dpi = 1200)
########

########
#Look at antigen specific populations d0-d15 for all 3 specificities
calcs <- naiveDF %>%
          group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
          summarize(n = n()) %>%
          filter(Timepoint %in% c("Day 0", "Day 15"))

calcs%>%
ggplot(aes(x = Timepoint, y = n))+
  geom_line(aes(group = Subject, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  scale_x_discrete(limits = c("Day 0", "Day 15"))+
  scale_color_manual(values = allColors3)+
  scale_fill_manual(values = allColors3)+
  facet_grid(cols = vars(Booster), rows = vars(adj.ProtoOmi), axes= "all_x")+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 6),
        legend.key.size = unit(0.8, "lines"))
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "AntigenPopulations_d1d15.png"), width = 4, height =5, dpi = 1200)

#do stats
stats <- naiveDF %>%
  group_by(Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  filter(Timepoint %in% c("Day 0", "Day 15")) %>%
  pivot_wider(names_from = Timepoint, values_from = n) %>% mutate(grouping = paste0(Booster, ", " ,adj.ProtoOmi)) %>%
  complete(adj.ProtoOmi, fill = list(`Day 0` = 0, `Day 15` = 0))

results <- by(stats, stats$grouping, function(x)
  wilcox.test(x$`Day 0`, x$`Day 15`, paired = TRUE))#[c(1:length(unique(calcs$Booster)))]
results <- type.convert(as.data.frame(do.call(rbind, results)), as.is=TRUE)
###############

#######
#Calculate SHM
ratpig <- citeRat %>%
            filter(!is.na(IgGQuant))

ratpigSHM <- observedMutations(ratpig, sequenceColumn = "sequence_alignment",
                                            germlineColumn = "germline_alignment",
                                            regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)

ratpigSHM2 <- ratpigSHM %>% mutate(DeNovo = BindingPopulation %in% c("BA.1-specific", "P-B+X+J+", "P-B+X+J-", "P-B+X-J+", "P-B+X-J-")) %>%
                filter(!(CorrectedBoost == "Prototype" & `Infection Status` == "Uninfected" & DeNovo))

ggplot(ratpigSHM2, aes(x= DeNovo, y= mu_freq))+
  geom_violin(aes(fill= DeNovo))+
  geom_boxplot(aes(fill = DeNovo), width = 0.2)+
  ylab("VH Mutation Frequency")+
  scale_fill_manual(values = c("TRUE" = "#C5A586", "FALSE" = "#EEE9E3"))+
  scale_x_discrete(limits = c("FALSE", "TRUE"), labels = c("Recall", "De Novo"))+
  theme_classic()+
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust= 1, vjust = 1),
        axis.title.x = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "Figure S5", "Denovo_Vs_Recall_SHM.png"), width = 2.5, height = 3)

denovo <- ratpigSHM2$mu_freq[ratpigSHM2$DeNovo]
recall <- ratpigSHM2$mu_freq[!ratpigSHM2$DeNovo]

wilcox.test(denovo, recall, paired= FALSE)
#######

#######
#write a table for sarah
naiveDF %>%
            group_by(InfectionRange, Booster, Subject, Timepoint, adj.ProtoOmi) %>%
            summarize(n = n()) %>% mutate(Proportion = n / sum(n)) %>% select(!n) %>%
            pivot_wider(names_from = Timepoint, values_from = Proportion) %>%
            mutate(InfectionRange = ifelse(is.na(InfectionRange), "Uninfected", InfectionRange))%>%
            select(InfectionRange, Booster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`) %>%
write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S5", "ProportionProbePositive_OverTime.xlsx"))

naiveDF %>%
  group_by(InfectionRange, Booster, Subject, Timepoint, adj.ProtoOmi) %>%
  summarize(n = n()) %>%
  pivot_wider(names_from = Timepoint, values_from = n) %>%
  mutate(InfectionRange = ifelse(is.na(InfectionRange), "Uninfected", InfectionRange))%>%
  select(InfectionRange, Booster, Subject, adj.ProtoOmi, `Day 0`, `Day 15`, `Day 90`, `Day 180`) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S5", "NumberProbePositive_OverTime.xlsx"))
#######