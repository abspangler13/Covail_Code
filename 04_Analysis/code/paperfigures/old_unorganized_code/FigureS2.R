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
#load in colors
allColors <- c("Prototype mRNA" = "#076894",
               "Prototype Protein" = "#73d4ff",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype/Beta mRNA" = "#c47002",
               "Prototype/Beta Protein" = "#f78c00",
               "Prototype/BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f",
               "Beta/BA.1 mRNA" = "wheat1")

#let's work on this
immunogenColors <- c("Prototype" = "#076894",
                     "Beta" = "#068041",
                     "Prototype/Beta" = "#c47002",
                     "Prototype/BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f",
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
         NotProto = rowSums(select(., c("Combined_Beta+/BA1+", "Combined_Beta+", "Combined_BA1+"))),
         Immunogen = str_replace_all(Immunogen, " \\+ ", "/"),
         Booster = str_replace_all(Booster, " \\+ ", "/"))

flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Prototype/Beta", "Beta", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1"))
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype Protein", 
                                                  "Prototype/Beta mRNA", "Prototype/Beta Protein", "Beta mRNA", "Beta Protein", 
                                                  "Prototype/BA.1 mRNA", "Omicron BA.1 mRNA", "Beta/BA.1 mRNA"))
flow$Timepoint <- factor(flow$Timepoint, levels = c("1", "15","90","180"))
#####

#####
#Compare between immunogens
####
#Beta
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype/Beta", "Beta"), Company != "Moderna") %>%
  group_by(Company, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoBeta / ProtoBeta [1],
         FoldProt = ProtoNotBeta / ProtoNotBeta[1]) %>%
  group_by(Company, Immunogen, Timepoint)%>%
  summarize(length = n(),
            meanC = mean(FoldCross),
            meanP = mean(FoldProt),
            seC = sd(FoldCross) / sqrt(length),
            seP = sd(FoldProt) / sqrt(length)
  )

#cross-reactive
ggplot(stats, aes(x = Timepoint, y=meanC, fill = Company))+
  geom_errorbar(aes(ymin = meanC-seC, ymax = meanC+seC, color = Company), width=0.2)+
  geom_line(aes(group = Company, color = Company), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Company))+
  ggtitle("Prototype+Beta+")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  #facet_grid(rows = vars(Platform))+
  ylim(0.5, 3)+
  geom_hline(yintercept = 1, linetype = 2)+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoBeta_byCompany.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoBeta_byCompany.svg"),width = 3, height = 1.8, device = "pdf")
dev.off()

#proto specific
ggplot(stats, aes(x = Timepoint, y=meanP, fill = Company))+
  geom_errorbar(aes(ymin = meanP-seP, ymax = meanP+seP, color = Company), width=0.2)+
  geom_line(aes(group = Company, color = Company), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Company))+
  ggtitle("Prototype+Beta-")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  #facet_grid(rows = vars(Platform))+
  ylim(0.5, 3)+
  geom_hline(yintercept = 1, linetype = 2)+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoOnly_byCompany.png"),width = 3, height = 1.8, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoOnly_byCompany.svg"),width = 3, height = 1.8, device = "pdf")
dev.off()

#write stats sheet
statistics <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype/Beta", "Beta"), Company != "Moderna") %>%
  group_by(Company, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoBeta / ProtoBeta [1],
         FoldProt = ProtoNotBeta / ProtoNotBeta[1]) %>%
  filter(Timepoint == 15) %>%
  select(Immunogen, Company, `Subject ID`, FoldCross, FoldProt)
write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "Figure S2", "ProtoBeta_CompanyFoldChangeComparison.csv"))


###########BA.1
stats <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1")) %>%
  group_by(Immunogen, Company, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoOmi / ProtoOmi[1],
         FoldProt = ProtoNotOmicron / ProtoNotOmicron[1],
         FoldBA1 = OmiNotProto / OmiNotProto[1]) %>%
  group_by(Immunogen, Company, Timepoint) %>%
  summarize(length = n(),
            meanC = mean(FoldCross),
            meanP = mean(FoldProt),
            seC = sd(FoldCross) / sqrt(length),
            seP = sd(FoldProt) / sqrt(length)
  )

#cross-reactive
ggplot(stats, aes(x = Timepoint, y=meanC, fill = Company))+
  geom_errorbar(aes(ymin = meanC-seC, ymax = meanC+seC, color = Company), width=0.2)+
  geom_line(aes(group = Company, color = Company), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Company))+
  ggtitle("Prototype+BA.1+")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  #facet_grid(rows = vars(Platform))+
  ylim(0.5, 5)+
  geom_hline(yintercept = 1, linetype = 2)+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoBA1_byCompany.png"),width = 4, height = 1.9, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoBA1_byCompany.svg"),width = 4, height = 1.9, device = "pdf")
dev.off()

#proto specific
ggplot(stats, aes(x = Timepoint, y=meanP, fill = Company))+
  geom_errorbar(aes(ymin = meanP-seP, ymax = meanP+seP, color = Company), width=0.2)+
  geom_line(aes(group = Company, color = Company), linewidth = 0.6)+
  geom_point(shape = 21, size =1.3, stroke = 0.6, aes(fill = Company))+
  ggtitle("Prototype+BA.1-")+
  ylab("Fold Change")+
  xlab("Timepoint")+
  scale_x_discrete(limits = c("1", "15", "90", "180"))+
  scale_fill_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  scale_color_manual(values = c("Moderna" = "darkblue", "Pfizer" = "lightblue", "Sanofi" = "#FCE036"))+
  #facet_grid(rows = vars(Platform))+
  ylim(0.5, 5)+
  geom_hline(yintercept = 1, linetype = 2)+
  facet_grid(cols = vars(Immunogen), labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(plot.title = element_text(size = 7, hjust = 0.5),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoOnly_notba1_byCompany.png"),width = 4, height = 1.9, units = "in", device = "png", dpi = 900)
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S2", "FoldCross_ProtoOnly_notba1_byCompany.svg"),width = 4, height = 1.9, device = "pdf")
dev.off()

#write stats sheet
statistics <- flow %>% filter(Immunogen %in% c("Prototype", "Prototype/BA.1", "Omicron BA.1", "Beta/BA.1")) %>%
  group_by(Company, Immunogen, `Subject ID`) %>%
  arrange(Timepoint) %>%
  mutate(FoldCross =  ProtoOmi / ProtoOmi [1],
         FoldProt = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
  filter(Timepoint == 15) %>%
  select(Immunogen, Company, `Subject ID`, FoldCross, FoldProt)
write.csv(statistics, here::here("04_Analysis", "data_objects", "paperfigures", "Figure S2", "ProtoBA1_CompanyFoldChangeComparison.csv"))

######
#reformatted boxplots
reformattedBA1 <- flow %>%
                    filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype/BA.1", "Beta/BA.1"), Timepoint %in% c("1", "15")) %>%
                    select(`Subject ID`, Timepoint, Immunogen, Platform, ProtoOmi, ProtoNotOmicron, OmiNotProto) %>%
                    pivot_longer(!c(`Subject ID`, Timepoint, Immunogen, Platform), names_to = "Population", values_to = "Total") %>%
                    mutate(Population = case_when(Population == "ProtoOmi" ~ "Cross-Reactive",
                                                  Population == "ProtoNotOmicron" ~ "Prototype-Specific",
                                                  Population == "OmiNotProto" ~ "BA.1-Specific"),
                           Population = factor(Population, levels = c("Cross-Reactive", "Prototype-Specific", "BA.1-Specific")),
                           Total = Total)

reformattedBA1 %>%
  ggplot(aes(x = Timepoint, y = Total))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.35, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("BA.1 Cross-Reactivity")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Population), rows = vars(Immunogen), axes = "all_x", labeller = label_wrap_gen(width = 10))+
  #ylim(c(-2.2,1.2))+
  scale_y_continuous(transform = "log10")+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BA1CrossReactivity_Wider.png"), width= 2.1, height =4, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BA1CrossReactivity_Wider.svg"), width= 2.1, height = 4)


#now for Beta
reformattedBeta <- flow %>%
  filter(Immunogen %in% c("Prototype", "Beta", "Prototype/Beta"), Timepoint %in% c("1", "15")) %>%
  select(`Subject ID`, Timepoint, Immunogen, Platform, ProtoBeta, ProtoNotBeta, BetaNotProto) %>%
  pivot_longer(!c(`Subject ID`, Timepoint, Immunogen, Platform), names_to = "Population", values_to = "Total") %>%
  mutate(Population = case_when(Population == "ProtoBeta" ~ "Cross-Reactive",
                                Population == "ProtoNotBeta" ~ "Prototype-Specific",
                                Population == "BetaNotProto" ~ "Beta-Specific"),
         Population = factor(Population, levels = c("Cross-Reactive", "Prototype-Specific", "Beta-Specific")),
         Total = Total)

reformattedBeta %>%
  ggplot(aes(x = Timepoint, y = Total))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.35, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Beta Cross-Reactivity")+
  #ylim(-2.2,1.2)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Population), rows = vars(Immunogen), axes = "all_x")+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BetaCrossReactivity_Wider.png"), width= 2.1, height =3.16, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BetaCrossReactivity_Wider.svg"), width= 2.1, height = 3.17)







########
#linear plots with full timecourse
#####
reformattedBA1 <- flow %>%
  filter(Immunogen %in% c("Prototype", "Omicron BA.1", "Prototype/BA.1", "Beta/BA.1")) %>%
  select(`Subject ID`, Timepoint, Immunogen, Platform, ProtoOmi, ProtoNotOmicron, OmiNotProto) %>%
  pivot_longer(!c(`Subject ID`, Timepoint, Immunogen, Platform), names_to = "Population", values_to = "Total") %>%
  mutate(Population = case_when(Population == "ProtoOmi" ~ "Cross-Reactive",
                                Population == "ProtoNotOmicron" ~ "Prototype-Specific",
                                Population == "OmiNotProto" ~ "BA.1-Specific"),
         Population = factor(Population, levels = c("Cross-Reactive", "Prototype-Specific", "BA.1-Specific")),
         Total = Total)

reformattedBA1 %>%
  ggplot(aes(x = Timepoint, y = Total))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("BA.1 Cross-Reactivity")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Population), rows = vars(Immunogen), axes = "all_x", labeller = label_wrap_gen(width = 10))+
  ylim(c(0,10))+
  #scale_y_continuous(transform = "log10")+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BA1CrossReactivity_Wider_linear.png"), width= 3, height =4, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BA1CrossReactivity_Wider_linear.svg"), width= 3, height = 4)


#now for Beta
reformattedBeta <- flow %>%
  filter(Immunogen %in% c("Prototype", "Beta", "Prototype/Beta")) %>%
  select(`Subject ID`, Timepoint, Immunogen, Platform, ProtoBeta, ProtoNotBeta, BetaNotProto) %>%
  pivot_longer(!c(`Subject ID`, Timepoint, Immunogen, Platform), names_to = "Population", values_to = "Total") %>%
  mutate(Population = case_when(Population == "ProtoBeta" ~ "Cross-Reactive",
                                Population == "ProtoNotBeta" ~ "Prototype-Specific",
                                Population == "BetaNotProto" ~ "Beta-Specific"),
         Population = factor(Population, levels = c("Cross-Reactive", "Prototype-Specific", "Beta-Specific")),
         Total = Total)

reformattedBeta %>%
  ggplot(aes(x = Timepoint, y = Total))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Beta Cross-Reactivity")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Population), rows = vars(Immunogen), axes = "all_x")+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BetaCrossReactivity_Wider_linear.png"), width= 3, height =3.16, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Prototype_BetaCrossReactivity_Wider_linear.svg"), width= 3, height = 3.17)

#################now do the non-wider boxplots
ba1 <- flow %>% filter(Immunogen %in% c("Omicron BA.1", "Prototype/BA.1", "Beta/BA.1", "Prototype"))
beta <- flow %>% filter(Immunogen %in% c("Beta", "Prototype", "Prototype/Beta"))

####make boxplots for BA1 populations
#proto omi
ba1 %>%
  ggplot(aes(x = Timepoint, y = ProtoOmi))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype+BA.1+")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "ProtoBA1_linear.png"), width= 3.2, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "ProtoBA1_linear.svg"), width= 3.2, height = 1.8)

#proto only
ba1 %>%
  ggplot(aes(x = Timepoint, y = ProtoNotOmicron))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype+BA.1-")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Proto_NotBA1_linear.png"), width= 3.2, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Proto_NotBA1_linear.svg"), width= 3.2, height = 1.8)

#BA.1 not prototype
ba1 %>%
  ggplot(aes(x = Timepoint, y = OmiNotProto))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype-BA.1+")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "BA1_notproto_linear.png"), width= 3.2, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "BA1_notproto_linear.svg"), width= 3.2, height = 1.8)

######it's beta time
#proto beta
beta %>%
  ggplot(aes(x = Timepoint, y = ProtoBeta))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype+Beta+")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "ProtoBeta_linear.png"), width= 2.5, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "ProtoBeta_linear.svg"), width= 2.5, height = 1.8)

#proto not beta
beta %>%
  ggplot(aes(x = Timepoint, y = ProtoNotBeta))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype+Beta-")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Proto_NotBeta_linear.png"), width= 2.5, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Proto_NotBeta_linear.svg"), width= 2.5, height = 1.8)

#beta not proto
beta %>%
  ggplot(aes(x = Timepoint, y = BetaNotProto))+
  geom_line(aes(group = `Subject ID`, color = Immunogen), lwd = 0.2, alpha = 0.6)+
  geom_boxplot(width= 0.5, lwd = 0.3, outlier.size = 0.9, outlier.shape = 21, outlier.stroke = 0.2, fatten = 1.4, aes(fill = Immunogen))+
  ylab("Percent of IgG+")+
  xlab("Days Post-Immunization")+
  ggtitle("Prototype-Beta+")+
  #ylim(-2.2,1.2)+
  ylim(0,10)+
  geom_hline(yintercept = 0.08, linetype = 2, linewidth = 0.6)+
  #scale_y_continuous(transform = "log10")+
  scale_fill_manual(values = immunogenColors)+
  scale_color_manual(values = immunogenColors)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(strip.background = element_blank(),
        text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 8),
        strip.text.y = element_text(size = 8, face = "bold"),
        legend.position = "none")
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Beta_NotProto_linear.png"), width= 2.5, height =1.8, dpi = 800)
ggsave(here::here('04_Analysis', "plots", "paperfigures", "Figure S2", "Beta_NotProto_linear.svg"), width= 2.5, height = 1.8)
