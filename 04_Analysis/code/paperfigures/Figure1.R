library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(cowplot)
library(RColorBrewer)

#set the paper's color schema
immunogenColors <- c("Prototype" = "#FBB042",
                     "Beta" = "#726658",
                     "Prototype/Beta" = "#BE1E2D",
                     "Prototype/BA.1" = "#1D75BC",
                     "Omicron BA.1" = "#2AB673",
                     "Beta/BA.1" = "wheat1")

#load in the flow data and make the result
flowRaw <- read_xlsx(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Filtered_COVAILDataset_250318_CORRECTED.xlsx"))

flow <- flowRaw %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         Immunogen = str_replace_all(Immunogen, " \\+ ", "/"),
         Booster = str_replace_all(Booster, " \\+ ", "/"))

flow$Timepoint <- factor(flow$Timepoint, levels = c("1", "15","90","180"))
#####

#####
#Generate numbers for intro graphic
####
stats <- flow %>%
  group_by(Immunogen, Platform, Timepoint) %>%
  summarize(n = length(unique(`Subject ID`)))
####

######
#Fig 1b: Total Response boxplot
flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1"))

ggplot(flow[flow$TotalRBD != 0,], aes(x = Timepoint, y=log10(TotalRBD)))+
  stat_boxplot(geom= 'errorbar', width = 0.2, lwd = 0.35)+
  geom_line(aes(group = `Subject ID`, color = Immunogen), alpha = 1, lwd = 0.1)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.25, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4)+
  ylab("log10 Total RBD+ Memory (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= immunogenColors)+
  scale_color_manual(values= immunogenColors)+
  facet_grid(cols = vars(Immunogen), rows = vars(Company), axes = "all", labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizervsModernavsSanofi_totalrbd.png"), width = 7, height = 3.8, units = "in", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizervsModernavsSanofi_totalrbd.svg"), width = 7, height = 3.8, units = "in")
dev.off()

#make stats sheet
statistics <-  flow %>% filter(infect_flag == "0" & Timepoint %in% c(1, 15)) %>%
  select(Company, Immunogen, `Subject ID`, TotalRBD, Timepoint) %>%
  mutate(TotalRBD = log10(TotalRBD)) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD)

write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "TotalRBD_EveryGroup_NAsRemoved_log.csv"))
#####

# #####
#Compare pfizer vs moderna vs sanofi as platforms for each immunogen
#make a delimiter for manufacturer, calculate the values to be graphed
stats <- flow %>%
  group_by(Immunogen, Company, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1],
         FoldTotalRBD = log2(FoldTotalRBD),
         Immunogen = factor(Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1", "Delta/BA.1"))) %>%
  group_by(Immunogen, Company, Timepoint) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            sd = sd(FoldTotalRBD) / sqrt(length) )

#make a dataframe summarizing the n for moderna vs pfizer
stats2 <- flow %>%
  group_by(Immunogen, Company, Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1],
         FoldTotalRBD = log2(FoldTotalRBD)) %>%
  #filter(FoldTotalRBD < 20) %>%
  mutate(Immunogen = factor(Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1", "Delta/BA.1"))) %>%
  summarize(n = n())

p1 <- ggplot(stats, aes(x = Timepoint, y=mean, fill = Company))+
  geom_hline(yintercept = 0, linetype = 2, lwd = 0.3)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Company), width=0.2, lwd=0.35)+
  geom_line(aes(group = Company, color = Company))+
  geom_point(shape = 21, aes(fill = Company), size = 1.3)+
  ggtitle("Total RBD: Pfizer vs Moderna")+
  ylab("log2 Fold Change")+
  xlab("Days Post-Immunization")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  geom_vline(xintercept = -Inf)+
  theme_classic()+
  #scale_y_continuous(breaks = c(0.5,1,3,5), limits = c(0.5,5))+
  theme(
    plot.title = element_blank(),
    axis.title.y = element_text(size=8),
    axis.title.x = element_text(size=8),
    axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face="bold"),
    panel.spacing = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.1, 'cm'),
    legend.title = element_text(size = 8, face = "bold"),
    legend.margin=margin(0,0,0,0),
    axis.line = element_line(linewidth=0.4),
    legend.position = "top")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))

#plot_grid(p1, p2, ncol = 1, align = "hv", rel_heights = c(1, 0.8))
p1
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizerVsModerna_lineplot_tabled_fold.png"), width = 7, height=2, units = "in", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizerVsModerna_lineplot_tabled_fold.svg"), width = 7, height=2)
dev.off()

#make a sheet for stats
stats <- flow %>%
  group_by(Immunogen, Company, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1],
         FoldTotalRBD = log2(FoldTotalRBD)) %>%
  filter(Timepoint == "15") %>%
  select(`Subject ID`, Timepoint, Immunogen, Company, FoldTotalRBD)

write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "ModernavsPfizervsSanofi_FC.csv"))
#####

#####
#instead of plotting fold change, do %igg difference
#make a delimiter for manufacturer
stats <- flow %>%
  group_by(Immunogen, Company, Timepoint) %>%
  mutate(Immunogen = factor(Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1", "Delta/BA.1"))) %>%
  summarize(length = n(),
            mean = mean(TotalRBD),
            sd = sd(TotalRBD) / sqrt(length))

#make a dataframe summarizing the n for moderna vs pfizer
stats2 <- flow %>%
  group_by(Immunogen, Company, Timepoint) %>%
  mutate(Immunogen = factor(Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1", "Delta/BA.1"))) %>%
  summarize(n = n())

p1 <- ggplot(stats, aes(x = Timepoint, y=mean, fill = Company))+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Company), width=0.2, lwd=0.35)+
  geom_line(aes(group = Company, color = Company))+
  geom_point(shape = 21, aes(fill = Company), size = 1.3)+
  ggtitle("Total RBD: Pfizer vs Moderna")+
  ylab("Total RBD (% IgG)")+
  xlab("Days Post-Immunization")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  geom_vline(xintercept = -Inf)+
  theme_classic()+
  #scale_y_continuous(breaks = c(0.5,1,3,5), limits = c(0.5,5))+
  theme(
    plot.title = element_blank(), 
    axis.title.y = element_text(size=8),
    axis.title.x = element_text(size=8),
    axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=8),
    strip.background = element_blank(),
    strip.text = element_text(size = 8, face="bold"),
    panel.spacing = unit(0.6, "lines"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.1, 'cm'),
    legend.title = element_text(size = 8, face = "bold"),
    legend.margin=margin(0,0,0,0),
    axis.line = element_line(linewidth=0.4),
    legend.position = "top")+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))

p1
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizerVsModerna_lineplot_tabled_percentigg.png"), width = 7, height=2, units = "in", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "PfizerVsModerna_lineplot_tabled_percentigg.svg"), width = 7, height=2)
dev.off()

#make a sheet for stats
stats <- flow %>%
  filter(Timepoint %in% c("1","15")) %>%
  select(`Subject ID`, Timepoint, Immunogen, Company, TotalRBD)

write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "ModernavsPfizervsSanofi_percentigg.csv"))
#####

#####
#Compare by immunogen
#####
# #Fig 1d
# stats <- flow %>% filter(Infection == "N") %>% group_by(Immunogen, Company, `Subject ID`) %>%
#   arrange(Timepoint) %>%
#   mutate(FoldTotalRBD = TotalRBD / TotalRBD[1],
#          FoldTotalRBD = log2(FoldTotalRBD)) %>%
#   group_by(Immunogen, Company, Timepoint) %>%
#   summarize(length = n(),
#             mean = mean(FoldTotalRBD),
#             sd = sd(FoldTotalRBD)
#   )
# 
# ggplot(stats[stats$Company == "Pfizer",], aes(x = Timepoint, y=mean, fill = Immunogen))+ 
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
#   geom_line(aes(group = Immunogen, color = Immunogen))+
#   geom_point(shape = 21, aes(fill = Immunogen))+
#   ggtitle("Total RBD - Pfizer")+
#   ylab("log2 Fold Change")+
#   xlab("Days Post-Immunization")+
#   scale_x_discrete(limits = c("1", "15", "90", "180"))+
#   # scale_fill_manual(values = myColors)+
#   # scale_color_manual(values = myColors)+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   ylim(-1,2)+
#   #scale_y_continuous(breaks = c(0.3, 1, 2.0, 3), limits = c(0.3, 3))+
#   theme_classic()+
#   geom_hline(yintercept = 0, linetype = 2)+
#   theme(plot.title = element_text(size=9, hjust = 0.5), 
#         axis.title.y = element_text(size=8),
#         axis.title.x = element_text(size=8),
#         axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         panel.spacing = unit(0.6, "lines"),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.1, 'cm'),
#         legend.title = element_text(size = 7),
#         legend.margin=margin(0,0,0,0))+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Pfizer.png"),width = 2.4, height = 1.8, units = "in", device = "png", dpi = 600)
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Pfizer.svg"),width = 2.4, height = 1.8)
# dev.off()
# 
# ggplot(stats[stats$Company == "Moderna",], aes(x = Timepoint, y=mean, fill = Immunogen))+
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
#   geom_line(aes(group = Immunogen, color = Immunogen))+
#   geom_point(shape = 21, aes(fill = Immunogen))+
#   ggtitle("Total RBD - Moderna")+
#   ylab("log2 Fold Change")+
#   xlab("Days Post-Immunization")+
#   scale_x_discrete(limits = c("1", "15", "90", "180"))+
#   # scale_fill_manual(values = myColors)+
#   # scale_color_manual(values = myColors)+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   theme_classic()+
#   ylim(-1,2)+
#   geom_hline(yintercept = 0, linetype = 2)+
#   theme(plot.title = element_text(size=9, hjust=0.5), 
#         axis.title.y = element_text(size=8),
#         axis.title.x = element_text(size=8),
#         axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         panel.spacing = unit(0.6, "lines"),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.1, 'cm'),
#         legend.title = element_text(size = 7),
#         legend.margin=margin(0,0,0,0))+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Moderna.png"),width = 2.4, height = 1.8, units = "in", device = "png", dpi = 600)
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Moderna.svg"),width = 2.4, height = 1.8, units = "in")
# dev.off()
# 
# ggplot(stats[stats$Company == "Sanofi",], aes(x = Timepoint, y=mean, fill = Immunogen))+ #let's not include infection yet- we'll make that point later in this figure!
#   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
#   geom_line(aes(group = Immunogen, color = Immunogen))+
#   geom_point(shape = 21, aes(fill = Immunogen))+
#   ggtitle("Total RBD - Sanofi")+
#   ylab("log2 Fold Change")+
#   xlab("Days Post-Immunization")+
#   scale_x_discrete(limits = c("1", "15", "90", "180"))+
#   # scale_fill_manual(values = myColors)+
#   # scale_color_manual(values = myColors)+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   theme_classic()+
#   ylim(-1,1.75)+
#   geom_hline(yintercept = 0, linetype = 2)+
#   theme(plot.title = element_text(size=9, hjust=0.5), 
#         axis.title.y = element_text(size=8),
#         axis.title.x = element_text(size=8),
#         axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=8),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 8),
#         panel.spacing = unit(0.6, "lines"),
#         legend.text = element_text(size = 6),
#         legend.key.size = unit(0.1, 'cm'),
#         legend.title = element_text(size = 7),
#         legend.margin=margin(0,0,0,0))+
#   guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Sanofi.png"),width = 2.4, height = 1.8, units = "in", device = "png", dpi = 600)
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_FoldChangeTotalRBD_Sanofi.svg"),width = 2.4, height = 1.8, units = "in")
# dev.off()
# 
# #write a file for stats
# stats <- flow %>% filter(Infection == "N" & !Platform == "Protein") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
#   arrange(Timepoint) %>%
#   mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
#   select(Booster, Immunogen, `Subject ID`, Timepoint, FoldTotalRBD) %>%
#   pivot_wider(names_from = Timepoint, values_from = FoldTotalRBD)
# 
# write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "FCTotalRBD_ImmunogenComparison.csv"))
# #####

#####
#Do log10 percent IgG for figure 1E
#Fig 1d
stats <- flow %>% filter(Infection == "N") %>% group_by(Immunogen, Company, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(TotalRBD = log10(TotalRBD)) %>%
  group_by(Immunogen, Company, Timepoint) %>%
  summarize(length = n(),
            mean = mean(TotalRBD),
            sd = sd(TotalRBD) / sqrt(length)
  )

ggplot(stats[stats$Company == "Pfizer",], aes(x = Timepoint, y=mean, fill = Immunogen))+ 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen))+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ggtitle("Total RBD - Pfizer")+
  ylab("Total RBD+ (log10 % IgG)")+
  xlab("Days Post-Immunization")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  # scale_fill_manual(values = myColors)+
  # scale_color_manual(values = myColors)+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  #scale_y_continuous(breaks = c(0.3, 1, 2.0, 3), limits = c(0.3, 3))+
  ylim(-1,1)+
  theme_classic()+
  theme(plot.title = element_text(size=9, hjust = 0.5), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.6, "lines"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 7),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_log10TotalRBD_Pfizer.svg"),width = 2.4, height = 1.8)
dev.off()

ggplot(stats[stats$Company == "Moderna",], aes(x = Timepoint, y=mean, fill = Immunogen))+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen))+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ggtitle("Total RBD - Moderna")+
  ylab("Total RBD+ (log10 % IgG)")+
  xlab("Days Post-Immunization")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  # scale_fill_manual(values = myColors)+
  # scale_color_manual(values = myColors)+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  theme_classic()+
  ylim(-1,1)+
  theme(plot.title = element_text(size=9, hjust=0.5), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.6, "lines"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 7),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_log10TotalRBD_Moderna.svg"),width = 2.4, height = 1.8, units = "in")
dev.off()

ggplot(stats[stats$Company == "Sanofi",], aes(x = Timepoint, y=mean, fill = Immunogen))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen))+
  geom_point(shape = 21, aes(fill = Immunogen))+
  ggtitle("Total RBD - Sanofi")+
  ylab("Total RBD+ (log10 % IgG+)")+
  xlab("Days Post-Immunization")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  # scale_fill_manual(values = myColors)+
  # scale_color_manual(values = myColors)+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  theme_classic()+
  ylim(-1,1)+
  theme(plot.title = element_text(size=9, hjust=0.5), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8),
        panel.spacing = unit(0.6, "lines"),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.1, 'cm'),
        legend.title = element_text(size = 7),
        legend.margin=margin(0,0,0,0))+
  guides(fill = guide_legend(title = "Vaccination"), color = guide_legend(title = "Vaccination"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1e_log10TotalRBD_Sanofi.svg"),width = 2.4, height = 1.8, units = "in")
dev.off()

#write a file for stats
stats <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(TotalRBD = log10(TotalRBD)) %>%
  select(Company, Immunogen, `Subject ID`, Timepoint, TotalRBD) %>%
  pivot_wider(names_from = Timepoint, values_from = TotalRBD)

write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 1", "log10TotalRBD_ImmunogenComparison.csv"))
#####