library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(RColorBrewer)

#####
#load in colors
allColors <- c("Prototype mRNA" = "#076894",
               "Prototype Protein" = "#73d4ff",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#f78c00",
               "Prototype + BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f",
               "Beta + BA.1 mRNA" = "wheat1")

#let's work on this
immunogenColors <- c("Prototype" = "#076894",
                     "Beta" = "#068041",
                     "Prototype + Beta" = "#c47002",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f",
                     "Beta + BA.1" = "wheat1")

#load in the flow data and make the necessary variables
flowRaw <- read_xlsx(here::here("04_Analysis", "data_objects", "figure_testing", "Infected_at_Baseline_flowdata.xlsx"))

flow <- flowRaw %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         NotProto = rowSums(select(., c("Combined_Beta+/BA1+", "Combined_Beta+", "Combined_BA1+"))))

flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Prototype + Beta", "Beta", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1"))
#####

########################Lineplots
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype + Beta", "Beta"), Company != "Moderna") %>%
  group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoBeta / ProtoBeta [1],
         FoldProt = ProtoNotBeta / ProtoNotBeta[1]) %>%
  group_by(Immunogen, Timepoint)%>%
  summarize(length = n(),
            meanC = mean(FoldCross),
            meanP = mean(FoldProt),
            seC = sd(FoldCross) / sqrt(length),
            seP = sd(FoldProt) / sqrt(length)
  )

#plot
#fold cross reactive
ggplot(stats, aes(x = Timepoint, y=meanC, fill = Immunogen))+
  geom_errorbar(aes(ymin = meanC-seC, ymax = meanC+seC, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Immunogen))+
  ggtitle("Prototype+/Beta+")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  #facet_grid(rows = vars(Platform))+
  ylim(0.5, 3.1)+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldCross_ProtoBeta_COMBINEDMRNAPROT.png"),width = 1.35, height = 1.8, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldCross_ProtoBeta_COMBINEDMRNAPROT.svg"),width = 1.35, height = 1.8, device = "svg")
dev.off()

#fold prototype specific
ggplot(stats, aes(x = Timepoint, y=meanP, fill = Immunogen))+
  geom_errorbar(aes(ymin = meanP-seP, ymax = meanP+seP, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke=0.6, aes(fill = Immunogen))+
  ggtitle("Prototype+/Beta-")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  ggtitle("Prototype+/Beta-")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  ylim(0.5, 3.1)+
  geom_hline(yintercept = 1, linetype = 2)+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  #facet_grid(rows = vars(Platform))+
  theme_classic()+
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldPrototype_ProtoNotBeta_COMBINEDMRNAPROT.png"),width = 2.65, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldPrototype_ProtoNotBeta_COMBINEDMRNAPROT.svg"),width = 2.65, height = 1.8, device = "pdf")
dev.off()

#Omicron cohorts
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1")) %>%
  group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoOmi / ProtoOmi[1],
         FoldProt = ProtoNotOmicron / ProtoNotOmicron[1],
         FoldBA1 = OmiNotProto / OmiNotProto[1]) %>%
  group_by(Immunogen, Timepoint) %>%
  summarize(length = n(),
            meanC = mean(FoldCross),
            meanP = mean(FoldProt),
            seC = sd(FoldCross) / sqrt(length),
            seP = sd(FoldProt) / sqrt(length)
  )

#cross
ggplot(stats, aes(x = Timepoint, y=meanC, fill = Immunogen))+
  geom_errorbar(aes(ymin = meanC-seC, ymax = meanC+seC, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Immunogen))+
  ggtitle("Prototype+/BA.1+")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  ylim(0.5, 3.1)+
  geom_hline(yintercept = 1, linetype = 2)+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2h_FoldCross_ProtoOmicron_COMBINEDMRNAPROT.png"),width = 1.35, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2h_FoldCross_ProtoOmicron_COMBINEDMRNAPROT.svg"),width = 1.35, height = 1.8, device = "svg")
dev.off()

#proto sp
ggplot(stats, aes(x = Timepoint, y=meanP, fill = Immunogen))+
  geom_errorbar(aes(ymin = meanP-seP, ymax = meanP+seP, color = Immunogen), width=0.2)+
  geom_line(aes(group = Immunogen, color = Immunogen), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke=0.6, aes(fill = Immunogen))+
  ggtitle("Prototype+/BA.1-")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  ylim(0.5, 3.1)+
  geom_hline(yintercept = 1, linetype = 2)+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  #facet_grid(rows = vars(Platform))+
  theme_classic()+
  theme(plot.title = element_text(size = 9, hjust = 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        axis.line = element_line(size = 0.3),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2h_FoldPrototype_ProtoNotOmi_COMBINEDMRNAPROT.png"),width = 2.65, height = 1.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2h_FoldPrototype_ProtoNotOmi_COMBINEDMRNAPROT.svg"),width = 2.65, height = 1.8, device = "svg")
dev.off()

#write statistics sheet
#ba.1
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1")) %>%
  group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoOmi / ProtoOmi[1],
         FoldProt = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
  select(Immunogen, `Subject ID`, Timepoint, FoldCross, FoldProt) %>%
  filter(Timepoint == "15")
write.csv(stats, here::here("04_Analysis", "data_objects", "figure_testing", "Figure 2", "FoldChange_Day15_CrossReactiveProtoSpecific_BA1.csv"))

#beta
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype + Beta", "Beta"), Company != "Moderna") %>%
  group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoBeta / ProtoBeta[1],
         FoldProt = ProtoNotBeta / ProtoNotBeta[1]) %>%
  select(Immunogen, `Subject ID`, Timepoint, FoldCross, FoldProt) %>%
  filter(Timepoint == "15")
write.csv(stats, here::here("04_Analysis", "data_objects", "figure_testing", "Figure 2", "FoldChange_Day15_CrossReactiveProtoSpecific_Beta.csv"))
#####


########################barplots
stats <- flow %>%
  filter(Immunogen %in% c("Prototype", "Beta", "Prototype + Beta"), Timepoint %in% c(1, 15), Company != "Moderna") %>%
  mutate(`Prototype-Specific` = ProtoNotBeta / (TotalRBD - NotProto),
         `Cross-Reactive` = ProtoBeta / (TotalRBD - NotProto)) %>%
  select(Immunogen, `Subject ID`, Timepoint, `Prototype-Specific`, `Cross-Reactive`) %>%
  pivot_longer(!c("Immunogen", "Subject ID", "Timepoint"), names_to = "Proportion") %>%
  group_by(Immunogen, Timepoint, Proportion) %>%
  summarize(mean = mean(value),
            n= n(),
            sd = sd(value)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = Proportion), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  #scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  #scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  ylab("Proportion Prototype+/Beta+")+
  xlab("Days Post-Immunization")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2_ProtoBetaCrossReactiveBarplot_COMBINEDMRNAPROT.png"),width = 3, height = 2, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2_ProtoBetaCrossReactiveBarplot_COMBINEDMRNAPROT.svg"),width = 3, height = 2)
dev.off()

#omicron
stats <- flow %>%
  filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype + BA.1", "Beta + BA.1"), Timepoint %in% c(1, 15)) %>%
  mutate(`Prototype-Specific` = ProtoNotOmicron / (TotalRBD - NotProto),
         `Cross-Reactive` = ProtoOmi / (TotalRBD - NotProto)) %>%
  select(Immunogen, `Subject ID`, Timepoint, `Prototype-Specific`, `Cross-Reactive`) %>%
  pivot_longer(!c("Immunogen", "Subject ID", "Timepoint"), names_to = "Proportion") %>%
  group_by(Immunogen, Timepoint, Proportion) %>%
  summarize(mean = mean(value),
            n= n(),
            sd = sd(value)) %>%
  mutate(se = sd / sqrt(n),
         cumusum = 1 - cumsum(mean),
         cumusum = ifelse(cumusum < 0, 0, cumusum))

ggplot(stats) +
  geom_bar(aes(x=Timepoint, y=mean, fill = Proportion), stat="identity", position="stack", color="black", linewidth = 0.3) +
  geom_errorbar(aes(x=Timepoint, ymin=ifelse(cumusum-se < 0, 0, cumusum-se), ymax= cumusum+se), width=0.2, alpha=0.9) +
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  #scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  #scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"))+
  ylab("Proportion Prototype+/BA.1+")+
  xlab("Days Post-Immunization")+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6), 
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 7),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2_ProtoOmiCrossReactiveBarplot_COMBINEDMRNAPROT.png"),width = 3.2, height = 2, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Figure2_ProtoOmiCrossReactiveBarplot_COMBINEDMRNAPROT.svg"),width = 3.2, height = 2)
dev.off()

##########
#fold change boxplots
#proto omi
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

stats <- flow %>% filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype + BA.1", "Beta + BA.1")) %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(Fold =  RatioCrossOmi / RatioCrossOmi [1])

ggplot(stats[stats$Timepoint == "15",], aes(x= Immunogen, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.2, linewidth = 0.5, aes(color = Immunogen))+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.shape = NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Immunogen), size = 0.8, stroke = 0.3)+
  ylab("Fold Change")+
  scale_fill_manual(values= immunogenColors)+
  scale_color_manual(values= immunogenColors)+
  scale_x_discrete(limits = c("Prototype","Prototype + BA.1", "Omicron BA.1", "Beta + BA.1"))+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 6), limits = c(0,7))+
  ggtitle("BA.1 Cross-Reactive : Prototype-Only")+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8, hjust= 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldChange_RatioCrossOmiToProtoSpecific.png"),width = 2.5, height = 2.8, dpi = 800, units = "in")
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldChange_RatioCrossOmiToProtoSpecific.svg"),width = 2.5, height = 2.8)
dev.off()

#proto beta
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta

stats <- flow %>% filter(Immunogen %in% c("Prototype", "Beta", "Prototype + Beta"), Company != "Moderna") %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(Fold =  RatioCrossBeta / RatioCrossBeta [1]) %>%
  filter(!is.infinite(Fold))

ggplot(stats[stats$Timepoint == "15",], aes(x= Immunogen, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.2, aes(color = Immunogen), position = position_dodge(width = 0.5))+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.2, outlier.shape = NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Immunogen), size = 0.8, stroke = 0.3)+
  ylab("Fold Change")+
  ggtitle("Beta Cross-Reactive : Prototype-Only")+
  scale_fill_manual(values= immunogenColors)+
  scale_color_manual(values= immunogenColors)+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 6), limits = c(0,7))+
  geom_hline(yintercept = 1, linetype=2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8, hjust = 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldChange_RatioCrossBetaToProtoSpecific_AllVaxCombined.png"),width = 2.2, height = 2.8, units = "in", device = "png", dpi = 600)
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "FoldChange_RatioCrossBetaToProtoSpecific_AllVaxCombined.svg"),width = 2.2, height = 2.8)
dev.off()

########write a stats sheet
#ba.1
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype + BA.1", "Beta + BA.1")) %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldBA1 =  RatioCrossOmi / RatioCrossOmi [1]) %>%
  filter(!is.infinite(FoldBA1), Timepoint == "15") %>%
  select(Immunogen, `Subject ID`, FoldBA1)
write.csv(stats, here::here("04_Analysis", "data_objects", "figure_testing", "Figure 2", "BA1_RatioCrossReactive_to_Prototype.csv"))

#beta
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Beta", "Prototype + Beta"), Company != "Moderna") %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldBeta =  RatioCrossBeta / RatioCrossBeta [1]) %>%
  filter(!is.infinite(FoldBeta), Timepoint == "15") %>%
  select(Immunogen, `Subject ID`, FoldBeta)
write.csv(stats, here::here("04_Analysis", "data_objects", "figure_testing", "Figure 2", "Beta_RatioCrossReactive_to_Prototype.csv"))




#######################################
#Make a comparison between infected and uninfected donors
flowFull <- read_xlsx(here::here("04_Analysis", "data_objects", "figure_testing", "Infected_at_Baseline_flowdata_uninfectedkept.xlsx"))

flowFull <- flowFull %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         NotProto = rowSums(select(., c("Combined_Beta+/BA1+", "Combined_Beta+", "Combined_BA1+")))) #%>% filter(Dataset == "Delta Panel")

flowFull$Immunogen <- factor(flowFull$Immunogen, levels = c("Prototype", "Prototype + Beta", "Beta", "Prototype + BA.1", "Omicron BA.1", "Beta + BA.1"))


#compare fold changes between infected and uninfected donors
flowFull$RatioCrossOmi <- flowFull$ProtoOmi / flowFull$ProtoNotOmicron

stats <- flowFull %>% filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype + BA.1", "Beta + BA.1")) %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(Fold =  RatioCrossOmi / RatioCrossOmi [1]) %>% filter(Timepoint == "15")

#compare between infection at baseline
ggplot(stats, aes(x= Immunogen, y = Fold))+
  geom_boxplot(aes(fill = infect_baseline), width= 0.8, lwd = 0.2, outlier.shape = NA)+
  geom_point(shape = 21, aes(fill = infect_baseline, group = infect_baseline),size = 0.4, stroke = 0.3, position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8))+
  ylab("Fold Change")+
  scale_fill_manual(values= c("Y" = "#d902ee", "N" = "#ffd79d"))+
  scale_x_discrete(limits = c("Prototype","Prototype + BA.1", "Omicron BA.1", "Beta + BA.1"))+
  ggtitle("BA.1 Cross-Reactive : Prototype-Only")+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 6), limits = c(0,7))+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8, hjust= 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        axis.line = element_line(size = 0.3),
        legend.title =  element_text(size = 8))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "BA1_ComparingRatio_infbaseline.png"),width = 5, height = 4)

######beta
flowFull$RatioCrossBeta <- flowFull$ProtoBeta / flowFull$ProtoNotBeta

stats <- flowFull %>% filter(Immunogen %in% c("Prototype", "Beta", "Prototype + Beta", "Beta + BA.1")) %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(Fold =  RatioCrossBeta / RatioCrossBeta [1]) %>%
  filter(!is.infinite(Fold))

#compare between infection at baseline
ggplot(stats, aes(x= Immunogen, y = Fold))+
  geom_boxplot(aes(fill = infect_baseline), width= 0.8, lwd = 0.2, outlier.shape = NA)+
  geom_point(shape = 21, aes(fill = infect_baseline, group = infect_baseline),size = 0.4, stroke = 0.3, position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8))+
  ylab("Fold Change")+
  scale_fill_manual(values= c("Y" = "#d902ee", "N" = "#ffd79d"))+
  scale_x_discrete(limits = c("Prototype","Prototype + Beta", "Beta", "Beta + BA.1"))+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 6), limits = c(0,7))+
  ggtitle("Beta Cross-Reactive : Prototype-Only")+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(plot.title = element_text(size = 8, hjust= 0.5),
        axis.title.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=8),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, face = "bold"),
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "figure_testing", "Figure 2", "Beta_ComparingRatio_infbaseline.png"),width = 4.2, height = 4)
