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

crossDF <- df %>% filter(adj.ProtoOmi == "Proto+Omi+")

#read in evolution data from Abby
evolving <- read.csv(here::here("01_raw-data", "Evo_dat_Timepoint_uniform.csv"))

########make activated dataset specifically
#identify clones that are activated at day 15 and follow their SHM over time
activatedClones <- unique(crossDF$clone_subject_id[crossDF$ClusterLabel %in% c("Acute Activated","Intermediate") & crossDF$Timepoint == "Day 15"])

#plot changes in SHM only among activated clones over time
activated <- crossDF %>% filter(clone_subject_id %in% activatedClones, Infection == "N") %>% mutate(mu_freq = mu_freq * 100)
#####

#set the colors
allColors <- c("Omicron BA.1 mRNA" = "#7C1D6f", 
               "Prototype + Omicron BA.1 mRNA" = "#DC3977",
               "Prototype mRNA" = "#045275")

immunogenColors <- c("BA.1" = "#7C1D6f", 
                     "BA.1 And Prototype" = "#DC3977",
                     "Prototype" = "#045275")
#####

#####
#look at bulk shm
ggplot(activated, aes(x = Timepoint, y = mu_freq, fill = OfficialBooster))+
  geom_violin()+
  geom_boxplot(width = 0.2, outlier.size = 0)+
  scale_fill_manual(values = allColors)+
  facet_grid(cols = vars(Booster))+
  ylim(0, 12.5)+
  theme_classic()+
  scale_x_discrete(limits=  c("Day 0", "Day 15", "Day 90", "Day 180"))+
  theme(text = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust =1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "activated_shm", "bulkSHM.png"), width = 4, height = 2)

activated %>% select(Booster, Timepoint, mu_freq) %>%
  write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "activated_shm", "bulkSHM.xlsx"))
#####

#####
#plot out lineages activated at day 15
perL <- activated %>%
          group_by(Booster, clone_subject_id) %>%
          mutate(TimepointC = ifelse(Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", Timepoint)) %>%
          filter(length(unique(TimepointC)) >= 3) %>%
          group_by(Booster, clone_subject_id, TimepointC) %>%
          summarize(mean = mean(mu_freq))

ggplot(perL, aes(x= TimepointC, y = mean, fill = Booster))+
  geom_line(aes(group = clone_subject_id, color = Booster), alpha = 0.6, linewidth = 0.4)+
  geom_point(shape = 21)+
  scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90/180"))+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Booster))+
  theme_classic()+
  theme(text = element_text(size = 7),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "activated_shm", "SHMPerLineage.png"), width = 5, height = 2.5)

#write statistics
check <- activated %>%
  group_by(Booster, clone_subject_id) %>%
  mutate(TimepointC = ifelse(Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", Timepoint)) %>%
  filter(length(unique(TimepointC)) >= 3) %>%
  group_by(Booster, clone_subject_id, TimepointC) %>%
  summarize(mean = mean(mu_freq)) %>%
  pivot_wider(names_from = TimepointC, values_from = mean) %>%
write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "activated_shm", "perlineage_shm.xlsx"))
#####

#####
#Evolution plots
aClones <- unique(df$clone_subject_id[df$ClusterLabel %in% c("Acute Activated","Intermediate") & df$Timepoint == "Day 15"])

#plot changes in SHM only among activated clones over time
activated2 <- df %>% filter(clone_subject_id %in% aClones, Infection == "N") %>% mutate(mu_freq = mu_freq * 100)

activatedEv <- activated2 %>%
               mutate(subject_clone_id = paste0("s", Subject,"_", clone_id, "_1"))

actEv <- evolving %>% filter(clone_id %in% activatedEv$subject_clone_id)

ggplot(actEv)+
  geom_jitter(aes(fill = adj.ProtoOmi, alpha = sig, shape = adj.ProtoOmi, x = Booster, y= slope), width = 0.2)+
  scale_shape_manual(values = c("Proto+Omi+" = 21, "Proto+Omi-" = 22, "Proto-Omi+"= 24))+
  scale_fill_manual(values = c("Proto+Omi+" = "#FFA630", "Proto+Omi-" =  "#4DA1A9", "Proto-Omi+" = "#D7E8BA"))+
  #scale_y_continuous(position = "right")+
  ylab("Slope")+
  scale_alpha_discrete(guide = "none")+
  ggtitle("Lineage Evolution")+
  #geom_text(data = summary, mapping = aes(label = n, x = Booster, y = 0.0005), size = 3.5)+
  theme_classic()+
  guides(shape = guide_legend(nrow = 2))+
  theme(text = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust =1, vjust =1),
        axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_blank())
ggsave(here::here("04_Analysis", "plots", "paperfigures", "activated_shm", "ClonalEvolution_uninfected_activatedbyvax.png"), width = 1.6, height =2.8)
ggsave(here::here("04_Analysis", "plots", "paperfigures", "activated_shm", "ClonalEvolution_uninfected_activatedbyvax.svg"), width = 1.6, height =2.8)
dev.off()
#####


# #########
# ########
# #do proto-specific analysis on per-lineage basis?
# ########
# ########
# proto <- df %>% filter(Infection == "N" & adj.ProtoOmi == "Proto+Omi-")
# 
# activatedProtoClones <- unique(proto$clone_subject_id[proto$ClusterLabel %in% c("Acute Activated","Intermediate") & proto$Timepoint == "Day 15"])
# 
# #plot changes in SHM only among activated clones over time
# activatedProto <- proto  %>% filter(clone_subject_id %in% activatedProtoClones, Infection == "N") %>% mutate(mu_freq = mu_freq * 100)
# 
# perL <- activatedProto %>%
#   group_by(Booster, clone_subject_id) %>%
#   mutate(TimepointC = ifelse(Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", Timepoint)) %>%
#   filter(length(unique(TimepointC)) >= 3) %>%
#   group_by(Booster, clone_subject_id, TimepointC) %>%
#   summarize(mean = mean(mu_freq))
# 
# ggplot(perL, aes(x= TimepointC, y = mean, fill = Booster))+
#   geom_line(aes(group = clone_subject_id, color = Booster), alpha = 0.6, linewidth = 0.4)+
#   geom_point(shape = 21)+
#   scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90/180"))+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   facet_grid(cols = vars(Booster))+
#   theme_classic()+
#   theme(text = element_text(size = 7),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         strip.background = element_blank(),
#         legend.position = "none")
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "activated_shm", "SHMPerLineage.png"), width = 5, height = 2.5)
# 
# #write statistics
# check <- activated %>%
#   group_by(Booster, clone_subject_id) %>%
#   mutate(TimepointC = ifelse(Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", Timepoint)) %>%
#   filter(length(unique(TimepointC)) >= 3) %>%
#   group_by(Booster, clone_subject_id, TimepointC) %>%
#   summarize(mean = mean(mu_freq)) %>%
#   pivot_wider(names_from = TimepointC, values_from = mean) %>%
#   write_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "activated_shm", "perlineage_shm.xlsx"))
# #####
# 
# #####
