library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)

#####
#set colors
allColors <- c("Prototype mRNA" = "#045275",
               "Prototype Protein" = "#0a86bf",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#f78c00",
               "Prototype + BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f",
               "Beta + BA.1 mRNA" = "firebrick",
               "Delta + BA.1 mRNA" = "goldenrod1")

immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#7CCBA2",
                     "Prototype + Beta" = "#c47002",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f",
                     "Beta + BA.1" = "firebrick",
                     "Delta + BA.1" = "goldenrod1")

#####load the old dataset
# old <- read_xlsx(here::here("01_raw-data", "FlowData","Stage23DeltaPaneltoincludeD0D15.xlsx"), sheet = "newdata")
# old <- old[!duplicated(old),] #remove duplicated entries
# 
# #fix existing variables
# old$Timepoint <- old$`Time point Guess`
# old$infect_flag <- as.character(old$infect_flag)
# old$Stage <- as.double(old$Stage)
# old$`Date run` <- as.double(old$`Date run`)
# old$Subject <- old$`Subject ID`
# 
# #let's calculate all of the standard things for this dataset
# old$probeset <- "Old"
# 
# #####
# #Add in really old flow data using delta as a probe with all the weird vaccination groups
# reallyOld <- read_xlsx(here::here("01_raw-data", "FlowData", "Stage1_Analysis_SFA_olddatamaybe.xlsx"), sheet = "Sheet5")
# reallyOldMetadata <- read.csv(here::here("01_raw-data", "FlowData", "bcell_unblinding.csv"))
# reallyOldMetadata2 <- read_xlsx(here::here("01_raw-data", "FlowData", "Stage1_Analysis_SFA_olddatamaybe.xlsx"), sheet = "count")
# 
# reallyOld <- reallyOld %>% mutate(Subject = reallyOldMetadata2$`Subject#`[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
#                     Timepoint = reallyOldMetadata$timepoint[match(`Sample:`, reallyOldMetadata$sn)],
#                     TimepointCheck = reallyOldMetadata2$Timepoint[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
#                     Treatment = reallyOldMetadata$treatment[match(`Sample:`, reallyOldMetadata$sn)],
#                     TreatmentCheck = reallyOldMetadata2$Group[match(`Sample:`, reallyOldMetadata2$`Sample:`)],
#                     infect_baseline = reallyOldMetadata$infect_baseline[match(`Sample:`, reallyOldMetadata$sn)],
#                     probeset = "Really Old") %>%
#             filter(!is.na(Subject)) %>%
#             mutate(Treatment = case_when(TreatmentCheck == "T2" ~ "1 Dose Omicron + Prototype (Moderna)",
#                                          TreatmentCheck == "T3" ~ "1 Dose Omicron (Moderna)",
#                                          TreatmentCheck == "T4" ~ "1 Dose Beta + Omicron (Moderna)",
#                                          TreatmentCheck == "T5" ~ "1 Dose Prototype (Moderna)",
#                                          TRUE ~ Treatment)) %>%
#             rename(`Beta+` = `Live/IgG/IgG++/Omi-Proto neg/Beta | Freq. of IgG++`,
#                   `Delta+` = `Live/IgG/IgG++/Omi-Proto neg/Delta | Freq. of IgG++`,
#                   `BA1+`=`Live/IgG/IgG++/Omi+ Prot-/Omi | Freq. of IgG++`,
#                   `Proto+`=`Live/IgG/IgG++/Proto+ Omi-/Proto | Freq. of IgG++`,
#                   `Proto+/Beta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Beta | Freq. of IgG++`,
#                   `Proto+/Delta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Delta | Freq. of IgG++`,
#                   `Proto+/Beta+/Delta+`=`Live/IgG/IgG++/Proto+ Omi-/Proto-Delta-Beta | Freq. of IgG++`,
#                   `Proto+/BA1+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi | Freq. of IgG++`,
#                   `Proto+/Beta+/BA1+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Beta | Freq. of IgG++`,
#                   `Proto+/Beta+/BA1+/Delta+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Beta-Delta | Freq. of IgG++`,
#                   `Proto+/BA1+/Delta+`=`Live/IgG/IgG++/Proto-Omicron/Proto-Omi-Delta | Freq. of IgG++`,
#                   `Beta+/Delta+` = `Live/IgG/IgG++/Omi-Proto neg/Delta-Beta | Freq. of IgG++`,
#                   `Beta+/BA1+/Delta+`=`Live/IgG/IgG++/Omi+ Prot-/Omi-Beta-Delta | Freq. of IgG++`,
#                   `Beta+/BA1+`=`Live/IgG/IgG++/Omi+ Prot-/Omi-Beta | Freq. of IgG++`,
#                   `BA1+/Delta+` = `Live/IgG/IgG++/Omi+ Prot-/Omi-Delta | Freq. of IgG++`,
#                   SN = `Sample:`)
# 
# #multiply by 100 to do percentage of IgG
# # check <- reallyOld %>% mutate(across(`Live/IgG/IgG++/Omi-Proto neg/Beta | Freq. of IgG++`:`Live/IgG/IgG++/Omi+ Prot-/Omi-Delta | Freq. of IgG++`, ~ .* 100))
# #commenting out above because flow data looks like it's already been converted to percentage of IGG
# df <- reallyOld[,colnames(reallyOld) %in% colnames(old)]
# df <- rbind(df, old[,colnames(old) %in% colnames(df)])
# 
# #there are a bunch of missing values for the date and the timepoint- let's get them fixed
# # df <- df %>%
# #       group_by(Subject) %>%
# #       mutate(Timepoint = case_when(n() == 1 ~ "Only One Sample Run",
# #                                    "Day 1" %in% unique(Timepoint) & is.na(Timepoint) ~ "Day 15",
# #                                    "Day 15" %in% unique(Timepoint) & is.na(Timepoint) ~ "Day 1",
# #                                    TRUE ~ Timepoint
# #                                    ))
# #commented out the above for now. There are some subjects who had 1 timepoint run repeatedly. I have no idea
# #what to do with these because they could be mislabelled and run at day 57. I'm going to just exclude them for now
# 
# #exclude the subjects with 1 or 3 runs
# df <- df %>%
#       group_by(Subject) %>%
#       filter(n() != 3 & n() != 1)
# 
# #now we can calculate different metrics
# df <- df %>%
#         mutate(TotalRBD = rowSums(across(`Beta+`:`BA1+/Delta+`)),
#                VariantNotProto = rowSums(select(across(`Beta+`:`BA1+/Delta+`), !starts_with("Proto+"))),
#                NonProtoSpecific = rowSums(select(across(`Beta+`:`BA1+/Delta+`), -c(`Proto+`))),
#                OfficialBooster = case_when(Treatment == "1 Dose Beta + Omicron (Moderna)" ~ "Beta + BA.1 mRNA",
#                                            Treatment == "1 Dose Delta + Omicron (Moderna)" ~ "Delta + BA.1 mRNA",
#                                            Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
#                                            Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
#                                            Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
#                                            Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
#                                            Treatment == "Beta (Sanofi)" ~ "Beta Protein",
#                                            Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
#                                            Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
#                                            Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
#                                            TRUE ~ "We have a problem"),
#                Immunogen = str_extract(OfficialBooster, ".+(?=( mRNA)|( Protein))"),
#                Platform = str_extract(OfficialBooster, "(mRNA)|(Protein)"))

df <- read.csv(here::here("01_raw-data","FlowData","BCU SFA_sentEMMES_comparing.csv")) %>%
        select(c("SN","TESTDATE","TESTTYPE","infect_baseline","treatment","timepoint","days_since_covid_exposure",
                 "dataset","pre_study_vax_type","Subject","usubjid","LYMPHP","BCELLSP","IgGP"),contains("..")) %>%
        mutate(TotalRBD = rowSums(across(`Proto..`:`DeltaBeta..1`)),
               OfficialBooster = case_when(Treatment == "1 Dose Beta + Omicron (Moderna)" ~ "Beta + BA.1 mRNA",
                                           Treatment == "1 Dose Delta + Omicron (Moderna)" ~ "Delta + BA.1 mRNA",
                                           Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                                           Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                                           Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                                           Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                                           Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                                           Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                                           Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                                           Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                                           TRUE ~ "We have a problem"),
               Immunogen = str_extract(OfficialBooster, ".+(?=( mRNA)|( Protein))"),
               Platform = str_extract(OfficialBooster, "(mRNA)|(Protein)"))

#write a csv for sarah
write.csv(df, here::here("04_Analysis", "data_objects", "paperfigures","AllCovailData_DeltaPanel.csv"))
#####


#########
#########
#what if we didn't subdivide by platform and pooled mrna and protein?
ggplot(df[df$infect_baseline == "N",], aes(x = Timepoint, y=TotalRBD))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_line(aes(group = Subject, color = Immunogen), alpha = 0.3, lwd = 0.4)+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ylab("Total RBD+ (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= immunogenColors)+
  scale_color_manual(values= immunogenColors)+
  facet_grid(cols = vars(Immunogen), axes = "all", labeller = label_wrap_gen(10))+
  scale_y_log10()+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_TotalRBD_CombinedmRNAAndProtein.png"),width = 7.4, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#write a csv
stats <- df %>% filter(infect_baseline == "N") %>%
  select(Immunogen, Subject, TotalRBD, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD)

write.csv(stats, here::here("04_Analysis","data_objects","paperfigures", "stats", "Figure 0", "Figure0_TotalRBDOverTime.csv"))
#####

#####
#Figure 0d: Calculate fold change in total RBD
stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) #%>%
  # group_by(OfficialBooster, Immunogen, Platform, Timepoint) %>%
  # summarize(length = n(),
  #           mean = mean(FoldTotalRBD),
  #           se = sd(FoldTotalRBD) / sqrt(length)
  # )

ggplot(stats[stats$Timepoint != "Day 1",], aes(x = Immunogen, y=FoldTotalRBD))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Immunogen), size = 1, stroke = 0.3)+
  ggtitle("Total RBD+")+
  ylab("Fold Change")+
  scale_fill_manual(values = immunogenColors)+
  scale_y_continuous(transform="log10")+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=7),
    axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_FoldChangeTotalRBD_combinedplatforms.png"), width = 2.9, height = 2.5, units = "in", device = "png", dpi = 600)
dev.off()

#write for stats
statistics <- stats %>%
  ungroup() %>%
  select(Immunogen, Subject, FoldTotalRBD, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)

write.csv(statistics, here::here("04_Analysis","data_objects","paperfigures", "stats", "Figure 0", "Figure0_TotalRBDOverTime_FoldChange.csv"))
#####

#####
#Figure 0e: fold change in variant-responsive population
#####
#Figure 0d: Calculate fold change in total RBD
# stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
#   arrange(Timepoint) %>%
#   mutate(FoldFully = `Proto+/Beta+/BA1+/Delta+` / `Proto+/Beta+/BA1+/Delta+`[1]) #%>%
# # group_by(OfficialBooster, Immunogen, Platform, Timepoint) %>%
# # summarize(length = n(),
# #           mean = mean(FoldTotalRBD),
# #           se = sd(FoldTotalRBD) / sqrt(length)
# # )
# 
# ggplot(stats[stats$Timepoint != "Day 1" & stats$FoldFully < 12,], aes(x = Platform, y=FoldFully, fill = OfficialBooster))+
#   stat_boxplot(geom= 'errorbar', width = 0.2)+
#   geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
#   ggtitle("Prototype+/Beta+/Delta+/BA.1+")+
#   ylab("Fold Change")+
#   xlab("Timepoint")+
#   scale_fill_manual(values = allColors)+
#   facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
#   theme_classic()+
#   geom_hline(yintercept = 1, linetype = 2)+
#   theme(
#     plot.title = element_text(size=8), 
#     axis.title.y = element_text(size=7),
#     axis.title.x = element_text(size=7),
#     axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#     axis.text.y = element_text(size=7),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 6.8),
#     panel.spacing = unit(0.6, "lines"),
#     # legend.text = element_text(size = 6),
#     # legend.key.size = unit(0.1, 'cm'),
#     # legend.title = element_text(size = 7),
#     # legend.margin=margin(0,0,0,0),
#     legend.position = "none")+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_FoldChange_FullyCrossReactive.png"), width = 4.8, height = 2.4, units = "in", device = "png", dpi = 600)
# dev.off()
# 
# #write for stats
# statistics <- stats %>%
#   ungroup() %>%
#   select(OfficialBooster, Subject, FoldTotalRBD, Timepoint) %>%
#   pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)
# 
# write_xlsx(statistics, here::here("04_Analysis","data_objects","paperfigures", "stats", "Figure 0", "Figure0_TotalRBDOverTime_FoldChange.xlsx"))
#####

#####
#Fig 0e: shifts in fully cross-reactive population over time
#the commented-out figure above is fine, but it's a common plot I've used for this paper.
#the prototype-specific population is weirdly low, so I'm gonna avoid showing that for now (do I have the wrong excel file??)
# stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
#   arrange(Timepoint) %>%
#   mutate(PropFull = `Proto+/Beta+/BA1+/Delta+` / TotalRBD,
#          deltaFully = PropFull - PropFull[1]) #%>%
# # group_by(OfficialBooster, Immunogen, Platform, Timepoint) %>%
# # summarize(length = n(),
# #           mean = mean(FoldTotalRBD),
# #           se = sd(FoldTotalRBD) / sqrt(length)
# # )
# 
# ggplot(stats[stats$Timepoint != "Day 1",], aes(x = Platform, y=deltaFully, fill = OfficialBooster))+
#   stat_boxplot(geom= 'errorbar', width = 0.2)+
#   geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
#   ggtitle("Prototype+/Beta+/Delta+/BA.1+")+
#   ylab("Change in Proportion")+
#   xlab("Timepoint")+
#   scale_fill_manual(values = allColors)+
#   facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
#   theme_classic()+
#   geom_hline(yintercept = 0, linetype = 2)+
#   theme(
#     plot.title = element_text(size=8),
#     axis.title.y = element_text(size=7),
#     axis.title.x = element_text(size=7),
#     axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#     axis.text.y = element_text(size=7),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 6.6),
#     panel.spacing = unit(0.6, "lines"),
#     # legend.text = element_text(size = 6),
#     # legend.key.size = unit(0.1, 'cm'),
#     # legend.title = element_text(size = 7),
#     # legend.margin=margin(0,0,0,0),
#     legend.position = "none")+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_Delta_FullyCrossReactive.png"), width = 4.9, height = 2.4, units = "in", device = "png", dpi = 600)
# dev.off()
#####

#####
#Fig 0f: shifts in fully cross-reactive population over time
#the commented-out figure above is fine, but it's a common plot I've used for this paper.
#the prototype-specific population is weirdly low, so I'm gonna avoid showing that for now (do I have the wrong excel file??)
# stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
#   arrange(Timepoint) %>%
#   mutate(ProtoOmi = `Proto+/BA1+`+`Proto+/Beta+/BA1+`+`Proto+/Beta+/BA1+/Delta+`+`Proto+/BA1+/Delta+`,
#          deltaFully = PropFull - PropFull[1]) #%>%
# # group_by(OfficialBooster, Immunogen, Platform, Timepoint) %>%
# # summarize(length = n(),
# #           mean = mean(FoldTotalRBD),
# #           se = sd(FoldTotalRBD) / sqrt(length)
# # )
# 
# ggplot(stats[stats$Timepoint != "Day 1",], aes(x = Platform, y=deltaFully, fill = OfficialBooster))+
#   stat_boxplot(geom= 'errorbar', width = 0.2)+
#   geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
#   ggtitle("Prototype+/Beta+/Delta+/BA.1+")+
#   ylab("Change in Proportion")+
#   xlab("Timepoint")+
#   scale_fill_manual(values = allColors)+
#   facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
#   theme_classic()+
#   geom_hline(yintercept = 0, linetype = 2)+
#   theme(
#     plot.title = element_text(size=8),
#     axis.title.y = element_text(size=7),
#     axis.title.x = element_text(size=7),
#     axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#     axis.text.y = element_text(size=7),
#     strip.background = element_blank(),
#     strip.text = element_text(size = 6.6),
#     panel.spacing = unit(0.6, "lines"),
#     # legend.text = element_text(size = 6),
#     # legend.key.size = unit(0.1, 'cm'),
#     # legend.title = element_text(size = 7),
#     # legend.margin=margin(0,0,0,0),
#     legend.position = "none")+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_Delta_FullyCrossReactive.png"), width = 4.9, height = 2.4, units = "in", device = "png", dpi = 600)
# dev.off()
#####

#####
#Fig 0f: fold change in anything that can bind either omicron, beta, or delta for relevant comparisons
df <- df %>% ungroup() %>%
      mutate(BA1Binding = rowSums(select(.,`Proto+/BA1+`, `Proto+/Beta+/BA1+`, `Proto+/Beta+/BA1+/Delta+`, `Proto+/BA1+/Delta+`, `BA1+`, `BA1+/Delta+`, `Beta+/BA1+`, `Beta+/BA1+/Delta+`)),
             DeltaBinding = rowSums(select(.,`Delta+`, `Proto+/Delta+`, `Proto+/Beta+/Delta+`, `Proto+/Beta+/BA1+/Delta+`, `Proto+/BA1+/Delta+`, `BA1+/Delta+`, `Beta+/Delta+`, `Beta+/BA1+/Delta+`)),
             BetaBinding = rowSums(select(.,`Beta+`, `Proto+/Beta+`, `Proto+/Beta+/Delta+`, `Proto+/Beta+/BA1+`, `Proto+/Beta+/BA1+/Delta+`, `Beta+/BA1+`, `Beta+/Delta+`, `Beta+/BA1+/Delta+`)))

stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalBA1 = BA1Binding / BA1Binding[1],
         FoldTotalDelta = DeltaBinding / DeltaBinding[1],
         FoldTotalBeta = BetaBinding / BetaBinding[1])

#relevant groups for BA1: Beta + BA.1, Delta + BA.1, Prototype + BA.1, Prototype (all mRNA)
stats %>% filter(OfficialBooster %in% c("Beta + BA.1 mRNA", "Delta + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA"), Timepoint != "Day 1", FoldTotalBA1 < 12) %>%
ggplot(aes(x = Immunogen, y=FoldTotalBA1, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Total BA.1+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype", "Omicron BA.1", "Prototype + BA.1", "Beta + BA.1", "Delta + BA.1"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=7),
    axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_FCInTotalBA1.png"), width = 2, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for delta: Delta + BA.1, Prototype, BA.1 by itself
stats %>% filter(OfficialBooster %in% c("Delta + BA.1 mRNA", "Prototype mRNA"), Timepoint != "Day 1", FoldTotalDelta < 12) %>%
  ggplot(aes(x = Immunogen, y=FoldTotalDelta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Total Delta+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype", "Delta + BA.1"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=7),
    axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FCInTotalDelta.png"), width = 1.2, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for beta: Proto, Proto+Beta, Beta+BA1, Beta
stats %>% filter(OfficialBooster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Beta + BA.1 mRNA", "Beta mRNA"), Timepoint != "Day 1", FoldTotalBeta < 12) %>%
  ggplot(aes(x = Immunogen, y=FoldTotalBeta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.size = 0.4)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Total Beta+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype", "Prototype + Beta", "Beta + BA.1", "Beta"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=7),
    axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FCInTotalBeta.png"), width = 2, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Fig 0: variant-specific populations over time
df <- df %>% ungroup() %>%
  mutate(BA1Only = rowSums(select(.,`BA1+`, `BA1+/Delta+`, `Beta+/BA1+`, `Beta+/BA1+/Delta+`)),
         DeltaOnly = rowSums(select(.,`Delta+`, `BA1+/Delta+`, `Beta+/Delta+`, `Beta+/BA1+/Delta+`)),
         BetaOnly = rowSums(select(.,`Beta+`, `Beta+/BA1+`, `Beta+/Delta+`, `Beta+/BA1+/Delta+`)))

stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
  arrange(Timepoint) %>%
  mutate(FoldBA1 = BA1Only / BA1Only[1],
         FoldDelta = DeltaOnly / DeltaOnly[1],
         FoldBeta = BetaOnly / BetaOnly[1])

stats %>% filter(OfficialBooster %in% c("Beta + BA.1 mRNA", "Delta + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & !is.infinite(FoldBA1) & FoldBA1 <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldBA1, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Prototype-/BA.1+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Beta + BA.1 mRNA", "Delta + BA.1 mRNA"))+
  ylim(c(0,6))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_FC_BA1Restricted.png"), width = 1.8, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for delta: Delta + BA.1, Prototype, BA.1 by itself
stats %>% filter(OfficialBooster %in% c("Delta + BA.1 mRNA", "Prototype mRNA") & !is.infinite(FoldDelta) & FoldDelta <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldDelta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Prototype-/Delta+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Delta + BA.1 mRNA"))+
  ylim(c(0,6))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FC_DeltaRestricted.png"), width = 1.2, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for beta: Proto, Proto+Beta, Beta+BA1, Beta
stats %>% filter(OfficialBooster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta + BA.1 mRNA", "Beta mRNA", "Beta Protein") & !is.infinite(FoldBeta) & FoldBeta <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldBeta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Prototype-/Beta+")+
  ylab("Fold Change")+
  xlab("")+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta + BA.1 mRNA", "Beta mRNA", "Beta Protein"))+
  ylim(c(0,6))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FC_BetaRestricted.png"), width = 2.1, height = 2.45, units = "in", device = "png", dpi = 600)
dev.off()
#####


#####
#Plot the ratios for this dataset
df <- df %>% ungroup() %>%
  mutate(RatioCrossBA1 = rowSums(select(.,`Proto+/BA1+`, `Proto+/Beta+/BA1+`, `Proto+/Beta+/BA1+/Delta+`, `Proto+/BA1+/Delta+`)) / rowSums(select(.,`Proto+`, `Proto+/Beta+`, `Proto+/Delta+`, `Proto+/Beta+/Delta+`)),
         RatioCrossDelta = rowSums(select(.,`Proto+/Delta+`, `Proto+/Beta+/Delta+`, `Proto+/Beta+/BA1+/Delta+`, `Proto+/BA1+/Delta+`)) / rowSums(select(.,`Proto+`, `Proto+/Beta+`, `Proto+/BA1+`, `Proto+/Beta+/BA1+`)),
         RatioCrossBeta = rowSums(select(.,`Proto+/Beta+`, `Proto+/Beta+/Delta+`, `Proto+/Beta+/BA1+`, `Proto+/Beta+/BA1+/Delta+`)) / rowSums(select(.,`Proto+`, `Proto+/BA1+`, `Proto+/Delta+`, `Proto+/BA1+/Delta+`)))

stats <- df %>% filter(infect_baseline == "N") %>% group_by(OfficialBooster, Immunogen, Platform, Subject) %>%
  arrange(Timepoint) %>%
  mutate(FoldBA1 = RatioCrossBA1 / RatioCrossBA1[1],
         FoldDelta = RatioCrossDelta / RatioCrossDelta[1],
         FoldBeta = RatioCrossBeta / RatioCrossBeta[1])

stats %>% filter(OfficialBooster %in% c("Beta + BA.1 mRNA", "Delta + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & FoldBA1 <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldBA1, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Proto+/BA1+ : Proto+/BA1-")+
  ylab("Fold Change")+
  xlab("")+
  ylim(c(0,6))+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Beta + BA.1 mRNA", "Delta + BA.1 mRNA"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0d_FC_BA1_CrossRatio.png"), width = 1.8, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for delta: Delta + BA.1, Prototype
stats %>% filter(OfficialBooster %in% c("Delta + BA.1 mRNA", "Prototype mRNA") & FoldDelta <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldDelta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Proto/Delta:Proto")+
  ylab("Fold Change")+
  xlab("")+
  ylim(c(0,6))+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Delta + BA.1 mRNA"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FC_Delta_CrossRatio.png"), width = 1.2, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#relevant groups for beta: Proto, Proto+Beta, Beta+BA1, Beta
stats %>% filter(OfficialBooster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta + BA.1 mRNA", "Beta mRNA", "Beta Protein") & FoldBeta <= 6, Timepoint != "Day 1") %>%
  ggplot(aes(x = OfficialBooster, y=FoldBeta, fill = OfficialBooster))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = OfficialBooster), width= 0.5, lwd = 0.2, outlier.shape=NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = OfficialBooster), size = 1, stroke = 0.3)+
  ggtitle("Proto+/Beta+ : Proto+/Beta-")+
  ylab("Fold Change")+
  xlab("")+
  ylim(c(0,6))+
  scale_fill_manual(values = allColors)+
  scale_x_discrete(limits = c("Prototype mRNA", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta + BA.1 mRNA", "Beta mRNA", "Beta Protein"))+
  #facet_grid(cols = vars(Immunogen), scales = "free_x", space = "free_x", labeller = label_wrap_gen(10))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=8), 
    axis.title.y = element_text(size=7),
    axis.title.x = element_text(size=6.5),
    axis.text.x = element_text(size=6.8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=7),
    strip.background = element_blank(),
    strip.text = element_text(size = 6.8),
    panel.spacing = unit(0.6, "lines"),
    # legend.text = element_text(size = 6),
    # legend.key.size = unit(0.1, 'cm'),
    # legend.title = element_text(size = 7),
    # legend.margin=margin(0,0,0,0),
    legend.position = "none")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FC_Beta_CrossRatio.png"), width = 2.1, height = 2.45, units = "in", device = "png", dpi = 600)
dev.off()

statsList <- list()

statsList[["BA1"]] <- stats %>%
  ungroup() %>%
  filter(OfficialBooster %in% c("Prototype mRNA" ,"Beta + BA.1 mRNA", "Delta + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA") & FoldBA1 <= 6) %>%
  select(OfficialBooster, Subject, FoldBA1, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldBA1)

statsList[["Delta"]] <- stats %>%
  ungroup() %>%
  filter(OfficialBooster %in% c("Prototype mRNA", "Delta + BA.1 mRNA") & FoldDelta <= 6) %>%
  select(OfficialBooster, Subject, FoldDelta, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldDelta)

statsList[["Beta"]] <- stats %>%
  ungroup() %>%
  filter(OfficialBooster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Prototype + Beta Protein","Beta + BA.1 mRNA", "Beta mRNA", "Beta Protein") & FoldBeta <= 6) %>%
  select(OfficialBooster, Subject, FoldBeta, Timepoint) %>%
  pivot_wider(names_from = Timepoint, values_from = FoldBeta)

openxlsx::write.xlsx(statsList, here::here("04_Analysis","data_objects","paperfigures", "stats", "Figure 0", "Figure0_RatioOfCrossReactiveToPrototypeSpecific_FoldChange.xlsx"))

#####




########
########
#New plots
########
########
#barplots showing day 15 probe positive populations
#make a pivoted table
stats <- df %>% 
  select(`Beta+`:Timepoint, Subject, Immunogen, TotalRBD) %>%
  mutate(across(c(`Beta+`:`BA1+/Delta+`), .fns = ~ . / TotalRBD)) %>%
  pivot_longer(!c(Subject, Timepoint, Immunogen), names_to = "Population", values_to = "Proportion") %>%
  group_by(Immunogen, Timepoint, Population) %>%
  mutate(Population = factor(Population, levels = c(
    "Proto+/Beta+/BA1+/Delta+",
    "Proto+/Beta+/Delta+", "Proto+/Beta+", "Proto+/Delta+", "Proto+/Beta+/BA1+", "Proto+/BA1+/Delta+", "Proto+/BA1+",
    "Beta+", "Delta+", "BA1+", "Beta+/BA1+/Delta+", "Beta+/BA1+", "Beta+/Delta+", "Delta+/BA1+",
    "Proto+"
  ))) %>%
  summarize(mean = mean(Proportion),
            n = n(),
            sd = sd(Proportion)) %>%
  filter(Population != "TotalRBD") %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum),
         TypePopulation = case_when(Population == "Proto+/Beta+/BA1+/Delta+" ~ "Extensively Cross-Reactive",
                                    Population %in% c("Proto+/BA1+", "Proto+/BA1+/Delta+", "Proto+/Beta+", "Proto+/Beta+/BA1+", "Proto+/Delta+", "Proto+/Beta+/Delta+") ~ "Cross-Reactive",
                                    Population == "Proto+" ~ "Prototype-Specific",
                                    TRUE ~ "Variant-Specific"),
         TypePopulation = factor(TypePopulation, levels = c("Extensively Cross-Reactive", "Cross-Reactive", "Variant-Specific", "Prototype-Specific")))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = TypePopulation), stat="identity", position="stack", color="black", linewidth = 0.1) +
  #geom_errorbar(aes(x=Timepoint, ymin=cumusum-se, ymax= ifelse(cumusum+se >= 1, 1, cumusum+se)), width=0.2, alpha=1) + #i can't get this to order properly when filling by TypePopulation
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  scale_x_discrete(limits = c("Day 1", "Day 15"))+
  scale_fill_manual(values = c("Extensively Cross-Reactive" = "#009392",
                                "Cross-Reactive" = "#9CCB86",
                                "Variant-Specific" = "#EEB479",
                               "Prototype-Specific" = "#CF597E"))+
  #scale_fill_brewer(palette = "Blues")+
  ylab("Proportion of Total RBD Response")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_ProbePositivePopulationsOverTime.png"),width = 4.2, height = 2.3, units = "in", device = "png", dpi = 600)
dev.off()

#######
#Fold change figure
stats <- df %>% 
  select(`Beta+`:Timepoint, Subject, Immunogen, TotalRBD) %>%
  mutate(across(c(`Beta+`:`BA1+/Delta+`), .fns = ~ . / TotalRBD)) %>%
  pivot_longer(!c(Subject, Timepoint, Immunogen), names_to = "Population", values_to = "Proportion") %>%
  mutate(Population = factor(Population, levels = c(
    "Proto+/Beta+/BA1+/Delta+",
    "Proto+/Beta+/Delta+", "Proto+/Beta+", "Proto+/Delta+", "Proto+/Beta+/BA1+", "Proto+/BA1+/Delta+", "Proto+/BA1+",
    "Beta+", "Delta+", "BA1+", "Beta+/BA1+/Delta+", "Beta+/BA1+", "Beta+/Delta+", "Delta+/BA1+",
    "Proto+"
  ))) %>%
  filter(Population != "TotalRBD") %>%
  arrange(Timepoint) %>%
  group_by(Immunogen, Subject, Population) %>%
  mutate(FoldChange = Proportion / Proportion[[1]]) %>%
  filter(!is.nan(FoldChange), !is.infinite(FoldChange), Timepoint != "Day 1") %>%
  group_by(Immunogen, Population)%>%
  summarize(median = median(FoldChange)) %>%
  mutate(TypePopulation = case_when(Population == "Proto+/Beta+/BA1+/Delta+" ~ "Extensively Cross-Reactive",
                                    Population %in% c("Proto+/BA1+", "Proto+/BA1+/Delta+", "Proto+/Beta+", "Proto+/Beta+/BA1+", "Proto+/Delta+", "Proto+/Beta+/Delta+") ~ "Cross-Reactive",
                                    Population == "Proto+" ~ "Prototype-Specific",
                                    TRUE ~ "Variant-Specific"),
         TypePopulation = factor(TypePopulation, levels = c("Extensively Cross-Reactive", "Cross-Reactive", "Variant-Specific", "Prototype-Specific")))

ggplot(stats, aes(x= Population, y= Immunogen))+
  geom_point(aes(size = abs(median - 1), fill = median - 1), shape = 21)+
  scale_fill_gradient2(midpoint = 0, high = "red", low = "blue")+
  scale_x_discrete(position = "top")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=0, vjust=0),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_ProbePositiveFoldChange.png"),width = 6, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()

#do regular fold change instead of proportional changes
stats <- df %>% 
  select(`Beta+`:Timepoint, Subject, Immunogen, TotalRBD) %>%
  #mutate(across(c(`Beta+`:`BA1+/Delta+`), .fns = ~ . / TotalRBD)) %>%
  pivot_longer(!c(Subject, Timepoint, Immunogen), names_to = "Population", values_to = "Proportion") %>%
  mutate(Population = factor(Population, levels = c(
    "Proto+/Beta+/BA1+/Delta+",
    "Proto+/Beta+/Delta+", "Proto+/Beta+", "Proto+/Delta+", "Proto+/Beta+/BA1+", "Proto+/BA1+/Delta+", "Proto+/BA1+",
    "Beta+", "Delta+", "BA1+", "Beta+/BA1+/Delta+", "Beta+/BA1+", "Beta+/Delta+", "Delta+/BA1+",
    "Proto+"
  ))) %>%
  filter(Population != "TotalRBD") %>%
  arrange(Timepoint) %>%
  group_by(Immunogen, Subject, Population) %>%
  mutate(FoldChange = Proportion / Proportion[[1]]) %>%
  filter(!is.nan(FoldChange), !is.infinite(FoldChange), Timepoint != "Day 1") %>%
  group_by(Immunogen, Population)%>%
  summarize(median = median(FoldChange)) %>%
  mutate(TypePopulation = case_when(Population == "Proto+/Beta+/BA1+/Delta+" ~ "Extensively Cross-Reactive",
                                    Population %in% c("Proto+/BA1+", "Proto+/BA1+/Delta+", "Proto+/Beta+", "Proto+/Beta+/BA1+", "Proto+/Delta+", "Proto+/Beta+/Delta+") ~ "Cross-Reactive",
                                    Population == "Proto+" ~ "Prototype-Specific",
                                    TRUE ~ "Variant-Specific"),
         TypePopulation = factor(TypePopulation, levels = c("Extensively Cross-Reactive", "Cross-Reactive", "Variant-Specific", "Prototype-Specific")))

ggplot(stats, aes(x= Population, y= Immunogen))+
  geom_point(aes(size = abs(median - 1), fill = median - 1), shape = 21)+
  scale_fill_gradient2(midpoint = 0, high = "red", low = "blue")+
  scale_x_discrete(position = "top")+
  ggtitle("Median Fold Changes in Binding Populations")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=8),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=0, vjust=0),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_ProbePositiveFoldChange_wholenumbers.png"),width = 6, height = 3, units = "in", device = "png", dpi = 600)

###evaluate specific pops
# stats <- df %>%
#   select(`Beta+`:Timepoint, Subject, Immunogen, TotalRBD) %>%
#   mutate(across(c(`Beta+`:`BA1+/Delta+`), .fns = ~ . / TotalRBD)) %>%
#   pivot_longer(!c(Subject, Timepoint, Immunogen), names_to = "Population", values_to = "Proportion") %>%
#   mutate(Population = factor(Population, levels = c(
#     "Proto+/Beta+/BA1+/Delta+",
#     "Proto+/Beta+/Delta+", "Proto+/Beta+", "Proto+/Delta+", "Proto+/Beta+/BA1+", "Proto+/BA1+/Delta+", "Proto+/BA1+",
#     "Beta+", "Delta+", "BA1+", "Beta+/BA1+/Delta+", "Beta+/BA1+", "Beta+/Delta+", "Delta+/BA1+",
#     "Proto+"
#   ))) %>%
#   filter(Population != "TotalRBD") %>%
#   arrange(Timepoint) %>%
#   group_by(Immunogen, Subject, Population) %>%
#   mutate(FoldChange = Proportion / Proportion[[1]]) %>%
#   filter(!is.nan(FoldChange), !is.infinite(FoldChange), Timepoint != "Day 1") %>%
#   ungroup()
# 
# ggplot(stats[stats$Population == "Proto+/Beta+/Delta+" & stats$FoldChange < 6,], aes(x = Immunogen, y= FoldChange))+
#   geom_point(aes(fill = Immunogen), shape = 21, position = position_jitterdodge())+
#   geom_hline(yintercept = 1)+
#   theme_classic()+
#   theme(legend.key.size = unit(0.2, 'cm'),
#         plot.title = element_text(size=6),
#         axis.title.y = element_text(size=7),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=7),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 6, face = "bold", margin = margin()),
#         panel.spacing = unit(0.3, "lines"),
#         legend.title = element_blank(),
#         legend.position = "bottom",
#         legend.text = element_text(size = 6),
#         legend.box.spacing = margin(0.5))

#######
#######
#Proto+ only and Extensively cross-reactive fold changes
#######
#######
stats <- df %>%
        group_by(Immunogen, Subject) %>%
        arrange(Timepoint) %>%
        select(OfficialBooster, Immunogen, Subject, Timepoint, `Proto+`, `Proto+/Beta+/BA1+/Delta+`) %>%
        mutate(FoldProtoOnly = `Proto+` / `Proto+`[[1]],
               FoldExtensive = `Proto+/Beta+/BA1+/Delta+` / `Proto+/Beta+/BA1+/Delta+`[[1]]) %>%
       filter(!is.nan(FoldProtoOnly), !is.infinite(FoldProtoOnly), Timepoint != "Day 1")
  

#proto only
ggplot(stats, aes(x = Immunogen, y = FoldProtoOnly, fill = Immunogen))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.size = 0)+
  geom_jitter(shape = 21, width=0.03, size = 1, stroke = 0.3)+
  scale_fill_manual(values = immunogenColors)+
  geom_hline(yintercept = 1, linetype= 2)+
  ylab("Fold Prototype-Specific MBCs")+
  theme_classic()+
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FCInProtoOnly.png"),width = 3, height = 2.6, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(stats, aes(x = Immunogen, y = FoldExtensive, fill = Immunogen))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.size = 0)+
  geom_jitter(shape = 21, width=0.03, size = 1, stroke = 0.3)+
  scale_fill_manual(values = immunogenColors)+
  geom_hline(yintercept = 1, linetype= 2)+
  ylab("Fold Broadly-Binding MBCs")+
  ylim(c(0,8))+
  theme_classic()+
  theme(legend.key.size = unit(0.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size=7),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 0", "Figure0_FCInExtensive.png"),width = 3, height = 2.6, units = "in", device = "png", dpi = 600)
dev.off()