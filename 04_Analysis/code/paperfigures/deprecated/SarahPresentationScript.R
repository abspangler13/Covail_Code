library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(rstatix)
library(gridExtra)
library(thematic)
library(readxl)
library(viridisLite)
library(stringdist)
library(corrplot)
library(ggpubr)
library(ggprism)
library(wesanderson)
library(writexl)
library(RColorBrewer)

#####
#set colors
#we'll be coloring by immunogen rather than platform now
# allColors <- c("Omicron BA.1 mRNA" = "#FF1F5B", 
#                  "Prototype + BA.1 mRNA" = "#009ADE",
#                  "Prototype mRNA" = "#FFC61E",
#                   "Prototype Protein" = "#ff6e54",
#                   "Prototype + Beta mRNA" = "#dd5182",
#                   "Prototype + Beta Protein" = "#955196",
#                   "Beta mRNA" = "#6573bf",
#                   "Beta Protein" = "#005e8a")

#based on Zissou1 color palette from wes anderson
allColors <- c("Prototype mRNA" = "#045275",
               "Prototype Protein" = "#0a86bf",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#0bd46c",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#fca532",
               "Prototype + BA.1 mRNA" = "#DC3977",
               "Omicron BA.1 mRNA" = "#7C1D6f")

immunogenColors <- c("Prototype" = "#045275",
                     "Beta" = "#7CCBA2",
                     "Prototype + Beta" = "#FCD39C",
                     "Prototype + BA.1" = "#DC3977",
                     "Omicron BA.1" = "#7C1D6f")

#load in the flow data
flow <- read_xlsx(here::here("01_raw-data", "FlowData","AllCOVAILMetadata_240314.xlsx"))
flow <- flow[!duplicated(flow),] #remove duplicated entries

#add in booster and remove anyone who has had an out-of-study boost
oosBoost <- unique(flow$`Subject ID`[flow$oosboost_flag == "Y"])
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

flow <- flow %>%
  filter(!`Subject ID` %in% oosBoost) %>%
  mutate(
    Treatment = ifelse(Treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", Treatment),
    Booster = case_when(Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                        Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA",
                        Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                        Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                        Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                        Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                        Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                        Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                        Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                        Treatment == "Omicron (Pfizer 1)" ~ "Omicron BA.1 mRNA",
                        Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + BA.1 mRNA",
                        TRUE ~ "Error"
    ),
    Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
  filter(`Subject ID` != 5752524948) #this donor has extraordinarily few memory B cells, so we should remove them

#set timepoint as character and then set order by declaring as factor
flow$`Time point Guess` <- as.character(flow$`Time point Guess`)
flow$`Time point Guess` <- factor(flow$`Time point Guess`, levels = c("1", "15", "90", "180"))

#set desired factor order for booster
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype Protein", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA",
                                                "Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein"))

#some funky stuff happened when Flavio exported from FlowJo, so we'll need to divide all of the values by 100 from the Pfizer groups
for(i in 1:nrow(flow)){
  if(flow$Treatment[i] %in% c("Omicron (Pfizer 1)", "Omicron + Wildtype/Prototype (Pfizer 1)")){
    for(j in c(seq(from = 24, to = 54, by =2))){
      flow[i,j] <- flow[i,j] / 100
    }
  }
} #wow i slayed that

#Prepare probe metrics
#Proto+Beta+
flow$`Proto+/Beta+/BA1+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+/BA1+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-BA1 | Freq. of IgG` * 100
flow$`Proto+/Beta+/XBB+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta-XBB | Freq. of IgG` * 100
flow$`Proto+/Beta+` <- flow$`Live/IgG/Proto+Beta+/Proto-Beta | Freq. of IgG` * 100

#Proto Only
flow$`Proto+` <- flow$`Live/IgG/Proto+/Proto | Freq. of IgG` * 100
flow$`Proto+/BA1+/XBB+` <- flow$`Live/IgG/Proto+/Proto-BA1-XBB | Freq. of IgG` * 100
flow$`Proto+/BA1+` <- flow$`Live/IgG/Proto+/Proto-BA1 | Freq. of IgG` * 100
flow$`Proto+/XBB+` <- flow$`Live/IgG/Proto+/Proto-XBB | Freq. of IgG` * 100

#Beta Only
flow$`Beta+` <- flow$`Live/IgG/Beta+/Beta | Freq. of IgG` * 100
flow$`Beta+/BA1+/XBB+` <- flow$`Live/IgG/Beta+/Beta-BA1-XBB | Freq. of IgG` * 100
flow$`Beta+/BA+` <- flow$`Live/IgG/Beta+/Beta-BA1 | Freq. of IgG` * 100
flow$`Beta+/XBB+` <- flow$`Live/IgG/Beta+/Beta-XBB | Freq. of IgG` * 100

#the m'crons
flow$`BA1+` <- flow$`Live/IgG/Proto neg Beta neg/BA1 | Freq. of IgG` * 100
flow$`XBB+` <- flow$`Live/IgG/Proto neg Beta neg/XBB | Freq. of IgG` * 100
flow$`BA1+/XBB+` <- flow$`Live/IgG/Proto neg Beta neg/BA1-XBB | Freq. of IgG` * 100

#Sum each of the totals:
#stage 1- we want a comparison between proto and omi-restricted compartments
flow$TotalRBD <- rowSums(flow[,c(57:71)])
flow$TotalBA1 <- rowSums(flow[,c(57,58,62,63,66,67,69,71)])
flow$ProtoNotBeta <- rowSums(flow[,c(61:64)])
flow$BetaNotProto <- rowSums(flow[,c(65:68)])
flow$ProtoBeta <- rowSums(flow[,c(57:60)])
flow$ProtoNotOmicron <- rowSums(flow[,c(59, 60, 61, 64)])
flow$OmicronNotPrototype <- rowSums(flow[,c(66, 67, 69, 71)])
flow$ProtoOmi <- rowSums(flow[,c(57, 58, 62, 63)])
flow$ProtoBetaNotOmi <- rowSums(flow[,c(59, 60)])
flow$ProtoOmiNotBeta <- rowSums(flow[,c(62, 63)])
flow$ProtoBetaOmiCrossReactive <- rowSums(flow[,c(57, 58)])
flow$ProtoBetaBA1NotXBB <- rowSums(flow[,c(58)])
flow$ProtoBetaBA1XBB <- rowSums(flow[,c(57)])
flow$ProtoBA1NotXBB <- rowSums(flow[,c(58, 63)])
flow$ProtoBA1XBB <- rowSums(flow[,c(57, 62)])

#set faceting variables now since doing so before would mess with column orders
flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster)

flow$infect_flag <- ifelse(is.na(flow$infect_flag), "0", flow$infect_flag)

#we will be merging our data with an earlier dataset run with a different probe panel (Beta on diff fluor, Delta variant instead of XBB)
#let's create a variable to differentiate the two
flow$probeset <- "New"

#####load the old dataset
old <- read_xlsx(here::here("01_raw-data", "FlowData","Stage23DeltaPaneltoincludeD0D15.xlsx"), sheet = "newdata")
old <- old[!duplicated(old),] #remove duplicated entries

#fix existing variables
old$`Time point Guess` <- factor(as.character(str_remove(old$`Time point Guess`, "Day ")), levels = c("1", "15", "90", "180"))
old$infect_flag <- as.character(old$infect_flag)
old$Stage <- as.double(old$Stage)
old$`Date run` <- as.double(old$`Date run`)

#let's calculate all of the standard things for this dataset
old$TotalRBD <- rowSums(old[,c(22:35)]) #I will exclude delta single-positives here as they are not relevant to any of our vaccination groups
old$TotalBA1 <- rowSums(old[,c(23, 24,26,27,30,31,34,35)])
old$ProtoNotBeta <- rowSums(old[,c(28,31, 32, 34)])
old$BetaNotProto <- rowSums(old[,c(22,24,25,27)])
old$ProtoBeta <- rowSums(old[,c(29,30,33,35)])
old$ProtoNotOmicron <- rowSums(old[,c(28,29,32,33)])
old$OmicronNotPrototype <- rowSums(old[,c(23,24,26,27)])
old$ProtoOmi <- rowSums(old[,c(30,31,34,35)])
old$ProtoBetaNotOmi <- rowSums(old[,c(29,33)])
old$ProtoOmiNotBeta <- rowSums(old[,c(31,34)])
old$ProtoBetaOmiCrossReactive <- rowSums(old[,c(35,30)])
old$probeset <- "Old"

#only keep columns present in flow
old <- old[,colnames(old) %in% colnames(flow)]

#merge with flow data
flow <- bind_rows(flow, old)

#rerun old variables so that the old data has them
infect <- unique(flow$`Subject ID`[flow$infect_flag == "Y" & !is.na(flow$infect_flag)])

flow <- flow %>%
  filter(!`Subject ID` %in% oosBoost) %>%
  mutate(
    Treatment = ifelse(Treatment == "1 Dose  Prototype (Moderna)", "1 Dose Prototype (Moderna)", Treatment),
    Booster = case_when(Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA",
                        Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA",
                        Treatment == "1 Dose Omicron (Moderna)" ~ "Omicron BA.1 mRNA",
                        Treatment == "1 Dose Omicron + Prototype (Moderna)" ~ "Prototype + BA.1 mRNA",
                        Treatment == "Beta (Pfizer 1)" ~ "Beta mRNA",
                        Treatment == "Beta + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + Beta mRNA",
                        Treatment == "Beta (Sanofi)" ~ "Beta Protein",
                        Treatment == "Beta + Prototype (Sanofi)" ~ "Prototype + Beta Protein",
                        Treatment == "Prototype (Sanofi)" ~ "Prototype Protein",
                        Treatment == "Omicron (Pfizer 1)" ~ "Omicron BA.1 mRNA",
                        Treatment == "Omicron + Wildtype/Prototype (Pfizer 1)" ~ "Prototype + BA.1 mRNA",
                        TRUE ~ "Error"
    ),
    Infection = ifelse(`Subject ID` %in% infect, "Y", "N"))%>%
  filter(`Subject ID` != 5752524948) #this donor has extraordinarily few memory B cells, so we should remove them

flow$Platform <- str_extract(flow$Booster, "(mRNA)|(Protein)")
flow$Immunogen <- case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster)

flow$Immunogen <- factor(flow$Immunogen, levels = c( "Omicron BA.1", "Prototype + BA.1", "Prototype", "Prototype + Beta", "Beta"))
#########

#####
#slide 5: total RBD response for just prototype vaccination
flow$TempDelim <- case_when(flow$Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA (Moderna)",
                            flow$Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA (Pfizer)",
                            flow$Treatment == "Prototype (Sanofi)" ~ "Prototype Protein (Sanofi)",
                            TRUE ~ "Irrelevant")

tempCols <- c(brewer.pal(5, "Blues")[5:3])

ggplot(flow[flow$infect_flag == "0" & flow$TotalRBD != 0 & flow$Booster %in% c("Prototype mRNA", "Prototype Protein"),], aes(x = `Time point Guess`, y=TotalRBD))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = TempDelim), alpha = 0.3)+
  geom_boxplot(aes(fill = TempDelim), width= 0.5)+
  ggtitle("Total RBD+ Memory B-Cell Counts Over Time")+
  ylab("Total RBD+ (Percentage of IgG+)")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values = tempCols)+
  scale_color_manual(values = tempCols)+
  facet_grid(cols = vars(TempDelim), labeller = label_wrap_gen(width = 15))+
  #ylim(c(0,13))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=13),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        panel.spacing = unit(0.4, "lines"))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide5_TotalRBD_ProtoVax.png"),width = 9, height = 4, units = "in", device = "png", dpi = 600)
dev.off()

#plot fold change together
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype Protein") & Infection == "N") %>% group_by(TempDelim, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(TempDelim, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

#plot
#Fold change Total Response- do we see a specific group dominating in total responses?
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = TempDelim))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = TempDelim), width=0.2)+
  geom_line(aes(group = TempDelim, color = TempDelim))+
  geom_point(shape = 21, aes(fill = TempDelim))+
  ggtitle("Mean Fold Change in Total RBD")+
  ylab("Mean Fold Change in Total RBD")+
  xlab("Timepoint")+
  scale_fill_manual(values = tempCols)+
  scale_color_manual(values = tempCols)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide5_MeanFC.png"),width = 7, height = 4, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#slide 6
flow$Booster <- factor(flow$Booster, levels = c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA", "Prototype Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta mRNA", "Beta Protein"))

ggplot(flow[flow$Booster %in% c("Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Prototype mRNA", "Prototype Protein") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=TotalRBD))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Total RBD+ Memory B-Cells Over Time by Flow")+
  ylab("Total RBD+ (Percentage of IgG+)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Platform), rows = vars(Immunogen), axes="all")+
  ylim(c(0,7))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide6_TotalRBDStage2.png"),width = 7, height =10, units = "in", device = "png", dpi = 600)
dev.off()

#fold change
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Prototype mRNA", "Prototype Protein") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Total RBD")+
  ylab("Mean Fold Change in Total RBD")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(rows = vars(Immunogen), axes = "all")+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide6_MeanFC.png"),width = 7, height = 10, units = "in", device = "png", dpi = 600)
dev.off()
######

######
#Slide 7
ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=TotalRBD))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Total RBD+ Memory B-Cells Over Time by Flow")+
  ylab("Total RBD+ (Percentage of IgG+)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,7))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide7_TotalRBD_AllmRNAImmunogens.png"),width = 14, height =5, units = "in", device = "png", dpi = 600)
dev.off()

#fold
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalRBD),
            se = sd(FoldTotalRBD) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Total RBD")+
  ylab("Mean Fold Change in Total RBD")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide7_MeanFC.png"),width = 7, height = 8, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Slide 8
ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=TotalBA1))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Total BA.1+ Memory B-Cells Over Time by Flow")+
  ylab("Total BA.1+ (Percentage of IgG+)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,7))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide8_TotalBA1_AllmRNAImmunogens.png"),width = 14, height =5, units = "in", device = "png", dpi = 600)
dev.off()

#fold
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalBA1 = TotalBA1 / TotalBA1[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldTotalBA1),
            se = sd(FoldTotalBA1) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Total BA1+ Memory B Cells")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide8_MeanFC_TotalBA1.png"),width = 7, height = 8, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Slide 9
uniqueVals <- unique(flow$Immunogen)
flow$Immunogen <- factor(flow$Immunogen, levels = c("Prototype", "Omicron BA.1", "Beta", "Prototype + BA.1","Prototype + Beta"))

#beta (stages 23)
ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Prototype-Beta Cross-Reactive Memory B-Cells Over Time")+
  ylab("Total Prototype+/Beta+ (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_protop_betap.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoNotBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Prototype-Positive, Beta-Negative Memory B-Cells")+
  ylab("Prototype+/Beta- (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_protop_betan.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=BetaNotProto))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Beta-Positive, Prototype-Negative Memory B-Cells")+
  ylab("Prototype-/Beta+ (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_proton_betap.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()

#stage1
ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoOmi))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Prototype-Omicron Cross-Reactive Memory B-Cells Over Time")+
  ylab("Total Prototype+/Omicron+ (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_protop_omip.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoNotOmicron))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Prototype-Positive, Omicron-Negative Memory B-Cells")+
  ylab("Prototype+/Omicron- (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_protop_omin.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=OmicronNotPrototype))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5)+
  ggtitle("Omicron-Positive, Prototype-Negative Memory B-Cells")+
  ylab("Prototype-/Omicron+ (Percentage of IgG)")+
  xlab("Timepoint")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Immunogen), axes="all", labeller = label_wrap_gen(width = 15))+
  ylim(c(0,6))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide9_proton_omip.png"),width = 10, height =5, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Slide 10
#fold change proto beta for stages 2/3
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoBeta = ProtoBeta / ProtoBeta[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldProtoBeta),
            se = sd(FoldProtoBeta) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Proto+/Beta+ Memory B Cells")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide10_MeanFC_ProtoBeta.png"),width = 9, height = 7, units = "in", device = "png", dpi = 600)
dev.off()

#now proto omi
stats <- flow %>% filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoOmi = ProtoOmi / ProtoOmi[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldProtoOmi),
            se = sd(FoldProtoOmi) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Proto+/Omicron+ Memory B Cells")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide10_MeanFC_ProtoOmicron.png"),width = 9, height = 7, units = "in", device = "png", dpi = 600)
dev.off()

#####

#####
#Slide 11
#do fold change for proto+ beta -
#fold change proto for stages 2/3
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoBeta = ProtoNotBeta / ProtoNotBeta[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldProtoBeta),
            se = sd(FoldProtoBeta) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Proto+/Beta- Memory B Cells")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide11_MeanFC_ProtoNotBeta.png"),width = 9, height = 5, units = "in", device = "png", dpi = 600)
dev.off()

#Do fold change from proto+omi-
stats <- flow %>% filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoOmi = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(FoldProtoOmi),
            se = sd(FoldProtoOmi) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster))+
  geom_point(shape = 21, aes(fill = Booster))+
  ggtitle("Mean Fold Change in Proto+/Omicron- Memory B Cells")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide11_MeanFC_ProtoNotOmicron.png"),width = 9, height = 5, units = "in", device = "png", dpi = 600)
dev.off()

###make the orig figure for proto-beta but exclude the protein groups
#Figure 1f. Relative ratios of cross-reactive to prototype-specific response
#calculate metrics
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

#remove donors with Day 0 or day 15 missing
cleanedMissing <- flow %>%
  group_by(`Subject ID`) %>%
  mutate(Missing = ifelse(length(intersect(`Time point Guess`, c("1", "15"))) < 2, "Yes", "No"),
         ContainsNARatio = ifelse(length(intersect(RatioCrossOmi, NA)) == 1 | length(intersect(RatioCrossBeta, NA)) == 1, "Yes","No"))%>%
  filter(Missing == "No")

#pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1D_BetaVaxVsProtoVax_CrossRatios.pdf"), width = 12, height = 4)
ggplot(cleanedMissing[cleanedMissing$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA",  "Prototype mRNA") & cleanedMissing$infect_flag == "0" & cleanedMissing$`Time point Guess` %in% c("1", "15") & cleanedMissing$ContainsNARatio == "No",], aes(x = `Time point Guess`, y=RatioCrossBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_line(aes(group = `Subject ID`))+
  geom_point(shape =21, aes(fill = Booster), size =2)+
  ggtitle("Ratio of Cross-Reactive to Prototype-Specific Response")+
  ylab("Prototype+/Beta+ / Prototype+Beta- Response")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(10))+
  ylim(c(0, 35))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=15), 
        axis.title.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "BetaVaxVsProtoVax_CrossRatios_proteinremoved.png"),width = 6, height = 5, units = "in", device = "png", dpi = 600)
dev.off()

#omi. now
#Omicron-specific
#pdf(file = here::here("04_Analysis", "plots", "paperfigures", "Figure 1", "Figure1D_OmicronVaxVsProtoVax_CrossRatios.pdf"), width = 9, height = 4)
ggplot(cleanedMissing[cleanedMissing$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & cleanedMissing$Infection == "N" & cleanedMissing$`Time point Guess` %in% c("1", "15") & cleanedMissing$ContainsNARatio == "No",], aes(x = `Time point Guess`, y=RatioCrossOmi))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_line(aes(group = `Subject ID`))+
  geom_point(shape =21, aes(fill = Booster), size=2)+
  ggtitle("Ratio Cross-Reactive to Prototype-Specific Response")+
  ylab("Ratio Proto+/Omicron+  /  Prototype+/Omicron- Response")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(10))+
  theme_classic()+
  theme(legend.key.size = unit(0.6, 'cm'),
        legend.text = element_text(size =12),
        plot.title = element_text(size=13), 
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=18,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=15),
        strip.background = element_blank(),
        strip.text = element_text(size = 15, face = "bold"),
        legend.position = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "OmicronVaxVsProtoVax_CrossRatios.png"),width = 6, height = 5, units = "in", device = "png", dpi = 600)
dev.off()

#####

#####
#Slide 12
#fold change in ratio of responses
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

stats <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldBetaRatio = RatioCrossBeta / RatioCrossBeta[1],
         FoldOmiRatio = RatioCrossOmi / RatioCrossOmi[1])

#plot
ggplot(stats[stats$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & stats$`Time point Guess` == "15",], aes(x = Immunogen, y=FoldOmiRatio, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_boxplot(aes(fill = Booster))+
  geom_point(shape =21, aes(fill = Booster))+
  ggtitle("FC in Cross-Reactive:Proto-Specific Ratio")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  ylim(c(0,5))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide12_MeanFC_ProtoOmiRatio.png"),width = 9, height = 5, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(stats[stats$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & stats$`Time point Guess` == "15",], aes(x = Immunogen, y=FoldBetaRatio, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_boxplot(aes(fill = Booster))+
  geom_point(shape =21, aes(fill = Booster))+
  ggtitle("FC in Cross-Reactive:Proto-Specific Ratio")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  ylim(c(0,5))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(
    plot.title = element_text(size=15), 
    axis.title.y = element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=14,angle = 45, hjust=1, vjust=1),
    axis.text.y = element_text(size=15),
    strip.background = element_blank(),
    strip.text = element_text(size = 15))+
  guides(fill = guide_legend(title="Vaccination"), color = "none")
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", "Slide12_MeanFC_ProtoBetaRatio.png"),width = 9, height = 5, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Statistics
#####

#####
#Slide 5
flow$TempDelim <- case_when(flow$Treatment == "1 Dose Prototype (Moderna)" ~ "Prototype mRNA (Moderna)",
                            flow$Treatment == "Wildtype/Prototype (Pfizer 1)" ~ "Prototype mRNA (Pfizer)",
                            flow$Treatment == "Prototype (Sanofi)" ~ "Prototype Protein (Sanofi)",
                            TRUE ~ "Irrelevant")

#let's do just day 0 to day 15 for now since this is the only part I can do paired that makes any sense
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype Protein") & Infection == "N") %>%
  select(`Subject ID`, TotalRBD, TempDelim, `Time point Guess`)    %>%     
  group_by(TempDelim, `Subject ID`) %>%
          pivot_wider(values_from = TotalRBD, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

#wilcox test
groups <- c("Proto Moderna", "Proto Pfizer", "Proto Protein")
results <- vector()
results[1] <- wilcox.test(stats$`1`[stats$TempDelim == "Prototype mRNA (Moderna)"],
                     stats$`15`[stats$TempDelim == "Prototype mRNA (Moderna)"], paired = TRUE)$p.value
results[2] <- wilcox.test(stats$`1`[stats$TempDelim == "Prototype mRNA (Pfizer)"],
                          stats$`15`[stats$TempDelim == "Prototype mRNA (Pfizer)"], paired = TRUE)$p.value
results[3] <- wilcox.test(stats$`1`[stats$TempDelim == "Prototype Protein (Sanofi)"],
                          stats$`15`[stats$TempDelim == "Prototype Protein (Sanofi)"], paired = TRUE)$p.value

print(data.frame(group = groups, result = results))

#FC
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype Protein") & Infection == "N") %>% group_by(TempDelim, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(TempDelim, `Time point Guess`) %>%
  select(`Subject ID`, FoldTotalRBD, TempDelim, `Time point Guess`) %>%
  pivot_wider(values_from = FoldTotalRBD, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`)) %>%
  ungroup()

stats %>% pairwise_wilcox_test( `15` ~ TempDelim, p.adjust.method = "bonferroni")
#####

#####Slide 6
#FC
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Beta Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Prototype mRNA", "Prototype Protein") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldTotalRBD, `Time point Guess`) %>%
  pivot_wider(values_from = FoldTotalRBD, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

results <- vector()
results[1] <- wilcox.test(stats$`15`[stats$Booster == "Prototype mRNA"],
                          stats$`15`[stats$Booster == "Prototype Protein"], paired = FALSE)$p.value
results[2] <- wilcox.test(stats$`15`[stats$Booster == "Prototype + Beta mRNA"],
                          stats$`15`[stats$Booster == "Prototype + Beta Protein"], paired = FALSE)$p.value
results[3] <- wilcox.test(stats$`180`[stats$Booster == "Beta Protein"],
                          stats$`180`[stats$Booster == "Beta mRNA"], paired = FALSE)$p.value
#####

#####Slide 7
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalRBD = TotalRBD / TotalRBD[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldTotalRBD, `Time point Guess`) %>%
  pivot_wider(values_from = FoldTotalRBD, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

write.csv(stats, here::here("04_Analysis", "data_objects", "paperfigures", "misc", "240328_SheetsForSarah","Slide7_FC.csv"))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                    wilcox.test(stats$`180`[stats$Booster == i],
                                stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
# day15 <- p.adjust(day15, method = "bonferroni")
# day90 <- p.adjust(day90, method = "bonferroni")
# day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)
####

#####Slide 8
#fold
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldTotalBA1 = TotalBA1 / TotalBA1[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldTotalBA1, `Time point Guess`) %>%
  pivot_wider(values_from = FoldTotalBA1, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                     wilcox.test(stats$`180`[stats$Booster == i],
                                 stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
day15 <- p.adjust(day15, method = "bonferroni")
day90 <- p.adjust(day90, method = "bonferroni")
day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)

#write data for Flavio
stats <- flow %>% filter(Infection == "N" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + Beta mRNA", "Prototype + BA.1 mRNA", "Beta mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, TotalBA1, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = TotalBA1) %>%
  na.omit()

write_xlsx(stats, here::here("04_Analysis","data_objects","paperfigures", "misc", "240328_SheetsForSarah", "Slide9_TotalBA1.xlsx"))


###


#####Slide 10
#prot beta stuff first
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoBeta = ProtoBeta / ProtoBeta[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldProtoBeta, `Time point Guess`) %>%
  pivot_wider(values_from = FoldProtoBeta, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                     wilcox.test(stats$`180`[stats$Booster == i],
                                 stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
# day15 <- p.adjust(day15, method = "bonferroni")
# day90 <- p.adjust(day90, method = "bonferroni")
# day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)

#now proto omi
stats <- flow %>% filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoOmi = ProtoOmi / ProtoOmi[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldProtoOmi, `Time point Guess`) %>%
  pivot_wider(values_from = FoldProtoOmi, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                     wilcox.test(stats$`180`[stats$Booster == i],
                                 stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
# day15 <- p.adjust(day15, method = "bonferroni")
# day90 <- p.adjust(day90, method = "bonferroni")
# day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)
###

#####Slide 11
stats <- flow %>% filter(Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoBeta = ProtoNotBeta / ProtoNotBeta[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldProtoBeta, `Time point Guess`) %>%
  pivot_wider(values_from = FoldProtoBeta, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                     wilcox.test(stats$`180`[stats$Booster == i],
                                 stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
# day15 <- p.adjust(day15, method = "bonferroni")
# day90 <- p.adjust(day90, method = "bonferroni")
# day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)

#Do fold change from proto+omi-
stats <- flow %>% filter(Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldProtoOmi = ProtoNotOmicron / ProtoNotOmicron[1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  select(`Subject ID`, FoldProtoOmi, `Time point Guess`) %>%
  pivot_wider(values_from = FoldProtoOmi, names_from = `Time point Guess`) %>%
  filter(!is.na(`1`) & !is.na(`15`))

unique <- unique(stats$Booster)
combo <- vector()
day15 <- vector()
day90 <- vector()
day180 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats$`15`[stats$Booster == i],
                                stats$`15`[stats$Booster == j], paired = FALSE)$p.value)
    
    day90 <- append(day90,
                    wilcox.test(stats$`90`[stats$Booster == i],
                                stats$`90`[stats$Booster == j], paired = FALSE)$p.value)
    
    day180 <- append(day180,
                     wilcox.test(stats$`180`[stats$Booster == i],
                                 stats$`180`[stats$Booster == j], paired = FALSE)$p.value)
  }
}

day15 <- p.adjust(day15, method = "bonferroni")
day90 <- p.adjust(day90, method = "bonferroni")
day180 <- p.adjust(day180, method = "bonferroni")

results <- data.frame(comparison = combo, day15 = day15, day90 = day90, day180 = day180)
####

#####Slide 12
#Slide 12
#fold change in ratio of responses
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

stats <- flow %>% filter(Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(FoldBetaRatio = RatioCrossBeta / RatioCrossBeta[1],
         FoldOmiRatio = RatioCrossOmi / RatioCrossOmi[1]) %>%
  filter(`Time point Guess` == "15")

stats23 <- stats %>% filter(Booster %in% c("Prototype mRNA", "Beta mRNA", "Prototype + Beta mRNA")) %>%
  select(`Subject ID`, FoldBetaRatio, `Time point Guess`) %>%
  pivot_wider(values_from = FoldBetaRatio, names_from = `Time point Guess`) %>%
  filter(!is.na(`15`)) %>%
  mutate(Booster = droplevels(Booster))

unique <- unique(stats23$Booster)
combo <- vector()
day15 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats23$`15`[stats$Booster == i],
                                stats23$`15`[stats$Booster == j], paired = FALSE)$p.value)
}
}
day15 <- p.adjust(day15, method = "bonferroni")
results <- data.frame(comparison = combo, day15 = day15)


stats1 <- stats %>% filter(Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>%
  select(`Subject ID`, FoldOmiRatio, `Time point Guess`) %>%
  pivot_wider(values_from = FoldOmiRatio, names_from = `Time point Guess`) %>%
  filter(!is.na(`15`)) %>%
  mutate(Booster = droplevels(Booster))

unique <- unique(stats1$Booster)
combo <- vector()
day15 <- vector()
for(i in unique){
  for(j in unique){
    combo <- append(combo, paste0(i, "_", j))
    
    day15 <- append(day15,
                    wilcox.test(stats1$`15`[stats$Booster == i],
                                stats1$`15`[stats$Booster == j], paired = FALSE)$p.value)
  }
}
day15 <- p.adjust(day15, method = "bonferroni")
results <- data.frame(comparison = combo, day15 = day15)

#####Slide 11
#make flavio a dataset with proto+variant- and proto-variant+
stats <- flow %>% filter(Infection == "N" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Beta mRNA", "Prototype + Beta mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, ProtoNotBeta, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = ProtoNotBeta) %>%
  na.omit()

write_xlsx(stats, here::here("04_Analysis","data_objects","paperfigures", "misc", "240328_SheetsForSarah", "ProtoBetaVax_PrototypePositive_BetaNegative.xlsx"))

stats <- flow %>% filter(Infection == "N" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Omicron BA.1 mRNA", "Prototype + BA.1 mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, OmicronNotPrototype, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = OmicronNotPrototype) %>%
  na.omit()

write_xlsx(stats, here::here("04_Analysis","data_objects","paperfigures", "misc", "240328_SheetsForSarah", "ProtoOmicronVax_PrototypeNegative_OmicronPositive.xlsx"))

stats <- flow %>% filter(Infection == "N" & `Time point Guess` %in% c("1", "15") & Booster %in% c("Prototype mRNA", "Beta mRNA", "Prototype + Beta mRNA")) %>%
  select(Treatment, `Subject ID`,probeset, BetaNotProto, `Time point Guess`) %>%
  pivot_wider(names_from = `Time point Guess`, values_from = BetaNotProto) %>%
  na.omit()

write_xlsx(stats, here::here("04_Analysis","data_objects","paperfigures", "misc", "240328_SheetsForSarah", "ProtoBetaVax_PrototypeNegative_BetaPositive.xlsx"))

#####




######
#CITESeq Data
#Slide 16/17 - clonal pie
#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_azimuth_ImmcantationRerunForPublicClones_demulti.rds"))
metadata <- seuObj@meta.data[seuObj@meta.data$Infection == "N",] %>% filter(ClusterLabel != "Naive")

splitDF <- split(metadata, metadata$Booster)
for(j in 1:length(splitDF)){
  
  filename <- paste0(names(splitDF)[j],"_NussenzweigStyleDonuts.pdf")
  df <- splitDF[[j]]
  
  pdf(file = here::here("04_Analysis", "plots", "paperfigures", "sarahpres", filename))
  for(i in unique(df$Subject)){
    placeholder <- df[df$Subject == i,]
    
    nonSinglets <- unique(placeholder$clone_subject_id[duplicated(placeholder$clone_subject_id) | duplicated(placeholder$clone_subject_id, fromLast=T)])
    placeholder$CloneStatus <- ifelse(placeholder$clone_subject_id %in% nonSinglets, placeholder$clone_subject_id, "Singlet")          
    
    placeholder$Timepoint <- ifelse(placeholder$Timepoint %in% c("Day 90", "Day 180"), "Day 90/180", placeholder$Timepoint)
    
    calcs <- placeholder %>%
      group_by(CloneStatus, Timepoint) %>%
      summarize(n = n()) %>%
      mutate(lab = case_when(CloneStatus == "Singlet" ~ "Singlet",
                             length(unique(Timepoint)) > 1 & "Day 0" %in% unique(Timepoint) ~ "Day 0 Expanded",
                             length(unique(Timepoint)) > 1 ~ "Expanded",
                             TRUE ~ "Single Timepoint"))
    
    placeholder <- placeholder %>%
      group_by(Timepoint, CloneStatus) %>%
      summarize(n= n()) %>%
      mutate(Proportion = n / sum(n),
             Total = sum(n),
             CloneStatus = fct_reorder(CloneStatus, Proportion, .desc=TRUE),
             adj.CloneStatus = case_when( CloneStatus == "Singlet" ~ "Singlet",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Day 0 Expanded"] ~ "Day 0 Expanded",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Expanded"] ~ "Expanded",
                                          CloneStatus %in% calcs$CloneStatus[calcs$lab == "Single Timepoint"] ~ "Single Timepoint"),
             adj.CloneStatus = fct(adj.CloneStatus, levels = c("Day 0 Expanded", "Expanded", "Single Timepoint", "Singlet")),
             Timepoint = factor(Timepoint, levels=c("Day 0", "Day 15", "Day 90/180")))%>%
      arrange(adj.CloneStatus)%>%
      mutate(ymax = cumsum(Proportion),
             ymin = c(0, head(ymax, n=-1)))
    
    #placeholder$CloneStatus <- factor(placeholder$CloneStatus, levels= c("Singlet","Expanded", "Day 0 Expanded", "Single Timepoint"))
    label <- c()
    label[1] <- unique(placeholder$Total[placeholder$Timepoint == "Day 0"])
    label[2] <- unique(placeholder$Total[placeholder$Timepoint == "Day 15"])
    label[3] <- unique(placeholder$Total[placeholder$Timepoint == "Day 90/180"])
    time <- c("Day 0", "Day 15", "Day 90/180")
    dat_text <- data.frame(label = label, Timepoint = time)
    dat_text$Timepoint <- factor(dat_text$Timepoint, levels = c("Day 0", "Day 15", "Day 90/180"))
    
    
    p <- ggplot(placeholder)+
      geom_rect(color= "black", linewidth=0.2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=adj.CloneStatus))+
      coord_polar(theta="y")+
      xlim(c(2,4))+
      scale_fill_manual(values = c("Singlet" = "#FFFFFF",
                                   "Day 0 Expanded" = "green4",
                                   "Expanded" = "#88CCEE",
                                   "Single Timepoint" = "gray80"))+
      ggtitle(paste0("Group, Subject: ",i, " ", names(splitDF)[j]))+
      guides(fill = "none")+
      facet_grid(cols=vars(Timepoint))+
      theme_void()+
      theme(plot.title = element_text(hjust=0.5))+
      geom_text(data = dat_text,
                mapping = aes(x=-Inf, y=-Inf, label = label),
                hjust = 0.5,
                vjust = 0.5,
                size = 10)
    
    print(p)
  }
  dev.off()
  rm(p)
}
rm(splitDF)
rm(p)
rm(calcs)
