library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)
library(cowplot)
library(RColorBrewer)

#####
#set colors
# allColors <- c("Prototype mRNA" = "#076894",
#                "Prototype Protein" = "#73d4ff",
#                "Beta mRNA" = "#068041",
#                "Beta Protein" = "#02ba5b",
#                "Prototype/Beta mRNA" = "#c47002",
#                "Prototype/Beta Protein" = "#f78c00",
#                "Prototype/BA.1 mRNA" = "#DC3977",
#                "Omicron BA.1 mRNA" = "#7C1D6f",
#                "Beta/BA.1 mRNA" = "wheat1")

#let's work on this
immunogenColors <- c("Prototype" = "#FBB042",
                     "Beta" = "#726658",
                     "Prototype/Beta" = "#BE1E2D",
                     "Prototype/BA.1" = "#1D75BC",
                     "BA.1" = "#2AB673",
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
         Booster = str_replace_all(Booster, " \\+ ", "/"),
         NotProto = rowSums(select(., c("Combined_Beta+/BA1+", "Combined_Beta+", "Combined_BA1+")))) %>%
  filter(Immunogen %in% c("Omicron BA.1", "Prototype/BA.1", "Prototype"), Platform == "mRNA") %>%
  mutate(Immunogen = str_remove(Immunogen, "Omicron "),
         Immunogen = factor(Immunogen, levels = c("Prototype", "Prototype/BA.1", "BA.1")))

flow$Timepoint <- factor(flow$Timepoint, levels = c("1", "15","90","180"))
#####

#####
#Make initial boxplot showing total anti-RBD response over time
##
ggplot(flow[flow$TotalRBD != 0,], aes(x = Timepoint, y=TotalRBD))+
  stat_boxplot(geom= 'errorbar', width = 0.2, lwd = 0.35)+
  geom_line(aes(group = `Subject ID`, color = Immunogen), alpha = 1, lwd = 0.1)+
  geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.25, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4)+
  ylab("Total RBD+ (% of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= immunogenColors)+
  scale_color_manual(values= immunogenColors)+
  scale_y_continuous(transform = "log10")+
  facet_grid(cols = vars(Immunogen), axes = "all_y")+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        axis.title.y = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.text.x = element_text(size=10,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=9),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        panel.spacing = unit(0.4, "lines"),
        legend.position = "none")
ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "TotalRBD_flow.svg"), width = 4.5, height = 2)
#####

#####
#make boxplots that show probe positive populations over time
stats <- flow %>% select(Immunogen, `Subject ID`, Timepoint, ProtoOmi, ProtoNotOmicron, OmiNotProto)%>%
  pivot_longer(cols = contains("Proto") ,names_to = "Population", values_to = "value") %>%
  mutate(Population = case_when(Population == "ProtoOmi" ~ "Prototype+BA.1+",
                                Population == "ProtoNotOmicron" ~ "Prototype+BA.1-",
                                Population == "OmiNotProto" ~ "Prototype-BA.1+",
                                TRUE ~ "Uh oh"),
         Population = factor(Population, levels = c("Prototype+BA.1+", "Prototype+BA.1-", "Prototype-BA.1+")))

#facet columns by population and rows by boost group to get the desired column
# ggplot(stats, aes(x = Timepoint, y = value))+
#   stat_boxplot(geom= 'errorbar', width = 0.2, lwd = 0.35)+
#   geom_line(aes(group = `Subject ID`, color = Immunogen), alpha = 1, lwd = 0.3)+
#   geom_boxplot(aes(fill = Immunogen), width= 0.5, lwd = 0.25, outlier.size = 1.2, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4)+
#   scale_fill_manual(values = immunogenColors)+
#   scale_color_manual(values = immunogenColors)+
#   ylab("% of IgG")+
#   facet_grid(rows = vars(Immunogen), cols = vars(Population), axes = "all_x", switch = "y")+
#   theme_classic()+
#   theme(text = element_text(size = 16),
#         legend.key.size = unit(0.6, 'cm'),
#         axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
#         strip.background = element_blank(),
#         panel.spacing = unit(0.4, "lines"),
#         legend.position = "none",
#         rect = element_rect(fill = "transparent"),
#         strip.placement = "outside",
#         strip.text = element_text(face = "bold"))
# ggsave(here::here("04_Analysis", "plots", "paperfigures", "keystone", "AntigenBindingPopulations_Flow.svg"), width = 6, height = 5.9)
#####

#####
#make barplots to show relative proportions of population
stats <- flow %>%
  filter(Immunogen %in% c("Prototype","BA.1", "Prototype/BA.1"), Timepoint %in% c(1, 15)) %>%
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
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(15), axes = "all_x")+
  #scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  #scale_x_discrete(limits = c("Day 0", "Day 15", "Day 90", "Day 180"))+
  scale_fill_manual(values = c("Prototype-Specific" = "lightgray",
                               "Cross-Reactive" = "#343148ff"),
                    labels = c("Prototype-Specific" = "Prototype+BA.1-",
                               "Cross-Reactive" = "Prototype+BA.1+"))+
  ylab("Proportion of Total RBD+")+
  xlab("Timepoint")+
  guides(fill = guide_legend(nrow = 2))+
  theme_classic()+
  theme(legend.key.size = unit(0.4, 'cm'),
        text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust=1, vjust=1),
        strip.background = element_blank(),
        # strip.text = element_text(size = 6.5, face = "bold", margin = margin()),
        # panel.spacing = unit(0.3, "lines"),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        legend.box.spacing = margin(0.5))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "FlowBarplots_CrossVsProto.svg"),width = 5.6, height = 3.2)
dev.off()
#####

#####
#Fold change ratio plot
##########
#fold change boxplots
#proto omi
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

stats <- flow %>% filter(Immunogen %in% c("Prototype", "BA.1", "Prototype/BA.1")) %>% group_by(Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(Fold =  RatioCrossOmi / RatioCrossOmi [1])

ggplot(stats[stats$Timepoint == "15",], aes(x= Immunogen, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.5, linewidth = 0.7)+
  geom_boxplot(aes(fill = Immunogen), width= 0.7, lwd = 0.5, outlier.shape = NA)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Immunogen), size = 1.3, stroke = 0.3)+
  ylab("Fold Change")+
  scale_fill_manual(values= immunogenColors)+
  scale_x_discrete(limits = c("Prototype","Prototype/BA.1", "BA.1"))+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 6), limits = c(0,7))+
  ggtitle("Ratio\nPrototype+BA.1+ / Prototype+BA.1-")+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(plot.title = element_text(size = 13, hjust= 0.5, color = "black"),
        axis.title.y = element_text(size=13, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13,angle = 45, hjust=1, vjust=1, color = "black"),
        axis.text.y = element_text(size=13, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold", color = "black"),
        legend.position = "none",
        axis.line = element_line(size = 0.4))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "keystone", "FoldChange_RatioCrossOmiToProtoSpecific.svg"),width = 3.2, height = 3.9)
dev.off()
#####