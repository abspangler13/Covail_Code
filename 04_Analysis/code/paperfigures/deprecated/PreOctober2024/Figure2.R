library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(gridExtra)
library(readxl)
library(stringr)
library(writexl)

#####
#set colors
allColors <- c("Prototype mRNA" = "#045275",
               "Prototype Protein" = "#0a86bf",
               "Beta mRNA" = "#068041",
               "Beta Protein" = "#02ba5b",
               "Prototype + Beta mRNA" = "#c47002",
               "Prototype + Beta Protein" = "#f78c00",
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
  filter(`Subject ID` != 5752524948, #this donor has extraordinarily few memory B cells, so we should remove them
         !(`Subject ID` == 5349564848 & `Time point Guess` == "15")) #this donor had problems with viability for day 15 sample)

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
flow$TotalRBD <- rowSums(flow[,c(57:69, 71)])
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
flow <- flow %>%
  mutate(Immunogen = case_when(flow$Booster %in% c("Prototype mRNA", "Prototype Protein") ~ "Prototype",
                            flow$Booster %in% c("Beta mRNA", "Beta Protein") ~ "Beta",
                            flow$Booster %in% c("Prototype + Beta mRNA", "Prototype + Beta Protein") ~ "Prototype + Beta",
                            flow$Booster == "Omicron BA.1 mRNA" ~ "Omicron BA.1",
                            flow$Booster == "Prototype + BA.1 mRNA" ~ "Prototype + BA.1",
                            TRUE ~ flow$Booster),
         Immunogen = factor(Immunogen, levels = c( "Omicron BA.1", "Prototype + BA.1", "Prototype", "Prototype + Beta", "Beta"))) %>%
        filter(`Subject ID` != 5450574848) #remove a specific donor for whom we don't have day 1 data for (but do have days 15 and 180)
#the donor removed above was not included in statistics
#####

#####
#Figure 2a. Showing differences in cross-reactivity towards different variants by booster
#Beta cross-reactive
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein"))

ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2, linewidth = 0.4)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Beta+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2a_BetaVaxVsProtoVax_ProtoBetaCrossReactive_mRNA.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Beta Protein", "Prototype + Beta Protein", "Prototype Protein") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2, linewidth = 0.4)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Beta+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2a_BetaVaxVsProtoVax_ProtoBetaCrossReactive_Protein.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

# ggplot(flow[flow$Booster %in% c("Beta Protein", "Prototype + Beta Protein", "Prototype Protein", "Prototype + Beta mRNA", "Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoBeta))+ #let's not include infection yet- we'll make that point later in this figure!
#   stat_boxplot(geom= 'errorbar', width = 0.2, linewidth = 0.4)+
#   geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
#   geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
#   ylab("Prototype+/Beta+ RBD")+
#   xlab("Days Post-Immunization")+
#   scale_fill_manual(values= allColors)+
#   scale_color_manual(values= allColors)+
#   facet_grid(cols = vars(Booster), labeller = label_wrap_gen(15))+
#   ylim(c(0,6))+
#   theme_classic()+
#   theme(axis.title.y = element_text(size=7),
#         axis.title.x = element_text(size=7),
#         axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
#         axis.text.y = element_text(size=7),
#         strip.background = element_blank(),
#         strip.text = element_text(size = 6, face = "bold"),
#         legend.position = "none",
#         axis.line = element_line(size = 0.3))
# ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2a_BetaVaxVsProtoVax_ProtoBetaCrossReactive.png"),width = 4.8, height = 2, units = "in", device = "png", dpi = 600)
# dev.off()

#Omicron-cross reactive
ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoOmi))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2, linewidth = 0.4)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Omicron+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2a_OmicronVaxVsProtoVax_ProtoOmiCrossReactive.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 2b. Prototype-specific RBD over time
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein"))

#beta group
ggplot(flow[flow$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoNotBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Beta- RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,1.5))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2b_BetaVaxVsProtoVax_ProtoSpec_mRNA.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Beta Protein",  "Prototype + Beta Protein", "Prototype Protein") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoNotBeta))+ #let's not include infection yet- we'll make that point later in this figure!
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Beta- RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,1.5))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2b_BetaVaxVsProtoVax_ProtoSpec_protein.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()


#Omicron group
ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=ProtoNotOmicron))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype+/Omicron- RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,5))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2b_OmicronVaxVsProtoVax_ProtoSpec.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 1g: Variant-specific (move to supp????)
flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta mRNA", "Beta Protein"))

#Beta-Specific
ggplot(flow[flow$Booster %in% c("Beta mRNA","Prototype + Beta mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=BetaNotProto))+ 
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype-/Beta+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2c_BetaVaxVsProtoVax_BetaSpec_mrna.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(flow[flow$Booster %in% c("Beta Protein", "Prototype + Beta Protein", "Prototype Protein") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=BetaNotProto))+ 
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype-/Beta+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2c_BetaVaxVsProtoVax_BetaSpec_protein.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

#Omicron spec
ggplot(flow[flow$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & flow$infect_flag == "0" & flow$TotalRBD != 0,], aes(x = `Time point Guess`, y=OmicronNotPrototype))+ 
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.3, lwd=0.3)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  ylab("Prototype-/Omicron+ RBD")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0,6))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2c_OmicronVaxVsProtoVax_OmicronSpec.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure 2. Relative ratios of cross-reactive to prototype-specific response
#calculate metrics
flow$RatioCrossBeta <- flow$ProtoBeta / flow$ProtoNotBeta
flow$RatioCrossOmi <- flow$ProtoOmi / flow$ProtoNotOmicron

flow$Booster <- factor(flow$Booster, levels = c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA", "Prototype Protein", "Prototype + Beta mRNA", "Prototype + Beta Protein", "Beta mRNA", "Beta Protein"))

#remove donors with Day 0 or day 15 missing
cleanedMissing <- flow %>%
  group_by(`Subject ID`) %>%
  mutate(Missing = ifelse(length(intersect(`Time point Guess`, c("1", "15"))) < 2, "Yes", "No"),
         ContainsNARatio = ifelse(length(intersect(RatioCrossOmi, NA)) == 1 | length(intersect(RatioCrossBeta, NA)) == 1, "Yes","No"))%>%
  filter(Missing == "No")

#beta
ggplot(cleanedMissing[cleanedMissing$Booster %in% c("Beta mRNA", "Prototype + Beta mRNA", "Prototype mRNA") & cleanedMissing$infect_flag == "0" & cleanedMissing$`Time point Guess` %in% c("1", "15") & cleanedMissing$ContainsNARatio == "No",], aes(x = `Time point Guess`, y=RatioCrossBeta))+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.6, lwd=0.3)+
  geom_point(shape =21, aes(fill = Booster), size =1.5, stroke = 0.3)+
  ylab("Prototype+Beta+ / Prototype+Beta-")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0, 35))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2d_BetaVaxVsProtoVax_CrossRatios_mrna.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(cleanedMissing[cleanedMissing$Booster %in% c("Beta Protein", "Prototype + Beta Protein", "Prototype Protein") & cleanedMissing$infect_flag == "0" & cleanedMissing$`Time point Guess` %in% c("1", "15") & cleanedMissing$ContainsNARatio == "No",], aes(x = `Time point Guess`, y=RatioCrossBeta))+
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.6, lwd=0.3)+
  geom_point(shape =21, aes(fill = Booster), size =1.5, stroke = 0.3)+
  ylab("Prototype+Beta+ / Prototype+Beta-")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  ylim(c(0, 35))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2d_BetaVaxVsProtoVax_CrossRatios_prot.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

#Omicron-specific
ggplot(cleanedMissing[cleanedMissing$Booster %in% c("Omicron BA.1 mRNA", "Prototype + BA.1 mRNA", "Prototype mRNA") & cleanedMissing$infect_flag == "0" & cleanedMissing$`Time point Guess` %in% c("1", "15") & cleanedMissing$ContainsNARatio == "No",], aes(x = `Time point Guess`, y=RatioCrossOmi))+ 
  geom_line(aes(group = `Subject ID`, color = Booster), alpha = 0.6, lwd=0.3)+
  geom_point(shape =21, aes(fill = Booster), size=1.5, stroke = 0.3)+
  ylab("Prototype+/BA.1+  /  Prototype+/BA.1-")+
  xlab("Days Post-Immunization")+
  scale_fill_manual(values= allColors)+
  scale_color_manual(values= allColors)+
  facet_grid(cols = vars(Booster), labeller = label_wrap_gen(13))+
  theme_classic()+
  theme(axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2d_OmicronVaxVsProtoVax_CrossRatios.png"),width = 2.4, height = 2, units = "in", device = "png", dpi = 600)
dev.off() #this turned out a lot better than I thought

#plot fold change just cuz
#proto omi
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & infect_flag == "0") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  RatioCrossOmi / RatioCrossOmi [1])

ggplot(stats[stats$`Time point Guess` == "15",], aes(x= Immunogen, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.2)+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Booster), size = 1, stroke = 0.3)+
  ylab("Fold Change")+
  scale_fill_manual(values= allColors)+
  scale_x_discrete(limits = c("Prototype","Prototype + BA.1", "Omicron BA.1"))+
  ylim(c(0.5,5.5))+
  ggtitle("Cross:Proto")+
  geom_hline(yintercept = 1, linetype = 2)+
  theme_classic()+
  theme(title = element_text(size = 6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "FoldChange_RatioCrossOmiToProtoSpecific.png"),width = 1.4, height = 2.2, units = "in", device = "png", dpi = 600)
dev.off()

#proto beta
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein") & infect_flag == "0") %>% group_by(Booster, Immunogen, Platform, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  RatioCrossBeta / RatioCrossBeta [1])

ggplot(stats[stats$`Time point Guess` == "15" & stats$Platform == "mRNA",], aes(x= Booster, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.2, aes(group = Booster), position = position_dodge(width = 0.5))+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Booster), size = 1, stroke = 0.3)+
  ylab("Fold Change")+
  ggtitle("Cross:Proto")+
  scale_fill_manual(values= allColors)+
  ylim(c(0,6))+ #an infected donor had a huge explosion in cross-reactive response at day 15, so I'm excluding them bc it seems really suspicious
  geom_hline(yintercept = 1, linetype=2)+
  theme_classic()+
  theme(title = element_text(size = 6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "FoldChange_RatioCrossBetaToProtoSpecific_mrna.png"),width = 1.4, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

ggplot(stats[stats$`Time point Guess` == "15" & stats$Platform == "Protein",], aes(x= Booster, y = Fold))+
  stat_boxplot(geom= 'errorbar', width = 0.2, aes(group = Booster), position = position_dodge(width = 0.5))+
  geom_boxplot(aes(fill = Booster), width= 0.5, lwd = 0.2, outlier.size = 0.3)+
  geom_jitter(shape = 21, width=0.025, aes(fill = Booster), size = 1, stroke = 0.3)+
  ylab("Fold Change")+
  ggtitle("Cross:Proto")+
  scale_fill_manual(values= allColors)+
  geom_hline(yintercept = 1, linetype = 2)+
  ylim(c(0,6))+ #an infected donor had a huge explosion in cross-reactive response at day 15, so I'm excluding them bc it seems really suspicious
  theme_classic()+
  theme(title = element_text(size = 6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "FoldChange_RatioCrossBetaToProtoSpecific_prot.png"),width = 1.4, height = 2.4, units = "in", device = "png", dpi = 600)
dev.off()

#####

#####
#Figure 1i: Fold changes for cross, proto-spec, and variant-specific(supp?) responses
#plot fold change together
#fold cross-reactive RBD
#stage 1
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & infect_flag == "0") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  ProtoOmi / ProtoOmi [1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(Fold),
            se = sd(Fold) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ 
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster), linewidth = 0.6)+
  geom_point(shape = 21, size =1.5, aes(fill = Booster))+
  ggtitle("Prototype+/BA.1+")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(title = element_text(size = 7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2e_FoldCross_ProtoOmi.png"),width = 1.8, height = 1.8, units = "in", device = "png", dpi = 600)
dev.off()

#stages 2/3
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  ProtoBeta / ProtoBeta [1]) %>%
  group_by(Booster, , Platform, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(Fold),
            se = sd(Fold) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster), linewidth = 0.6)+
  geom_point(shape = 21, size =1.5, aes(fill = Booster))+
  ggtitle("Prototype+/Beta+")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(rows = vars(Platform))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(title = element_text(size = 7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2h_FoldCross_ProtoBeta.png"),width = 1.8, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()
#####

#####
#Figure j: fold change for proto-specific population
#stage 1
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + BA.1 mRNA", "Omicron BA.1 mRNA") & infect_flag == "0") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  ProtoNotOmicron / ProtoNotOmicron [1]) %>%
  group_by(Booster, Immunogen, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(Fold),
            se = sd(Fold) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+ #let's not include infection yet- we'll make that point later in this figure!
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster), linewidth = 0.6)+
  geom_point(shape = 21, size =1.5, aes(fill = Booster))+
  ggtitle("Prototype+/BA.1-")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(title = element_text(size = 7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2f_FoldProtoSp_ProtoNotOmi.png"),width = 1.8, height = 1.8, units = "in", device = "png", dpi = 600)
dev.off()

#stages 2/3
stats <- flow %>% filter(Booster %in% c("Prototype mRNA", "Prototype + Beta mRNA", "Beta mRNA", "Prototype Protein", "Prototype + Beta Protein", "Beta Protein") & Infection == "N") %>% group_by(Booster, Immunogen, `Subject ID`) %>%
  arrange(`Time point Guess`) %>%
  mutate(Fold =  ProtoNotBeta / ProtoNotBeta [1]) %>%
  group_by(Booster, Immunogen, Platform, `Time point Guess`) %>%
  summarize(length = n(),
            mean = mean(Fold),
            se = sd(Fold) / sqrt(length)
  )

#plot
ggplot(stats, aes(x = `Time point Guess`, y=mean, fill = Booster))+
  geom_errorbar(aes(ymin = mean-se, ymax = mean+se, color = Booster), width=0.2)+
  geom_line(aes(group = Booster, color = Booster), linewidth = 0.6)+
  geom_point(shape = 21, size =1.5, aes(fill = Booster))+
  ggtitle("Prototype+/Beta-")+
  ylab("Mean Fold Change")+
  xlab("Timepoint")+
  scale_fill_manual(values = allColors)+
  scale_color_manual(values = allColors)+
  facet_grid(rows = vars(Platform))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype = 2)+
  theme(title = element_text(size = 7),
        axis.title.y = element_text(size=7),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7,angle = 45, hjust=1, vjust=1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold"),
        legend.position = "none",
        axis.line = element_line(size = 0.3))
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2f_FoldProtoSp_ProtoNotBeta.png"),width = 1.8, height = 2.8, units = "in", device = "png", dpi = 600)
dev.off()
#####

######
######
#extra figures
######
######
#showing what probe combinations exist among BA.1+ populations
stats <- flow %>% filter(probeset == "New") %>%
          select(`Subject ID`, Immunogen, `Time point Guess`, `BA1+`, `BA1+/XBB+`, 
                 `Proto+/Beta+/BA1+/XBB+`, `Proto+/Beta+/BA1+`, `Proto+/BA1+/XBB+`, `Proto+/BA1+`,
                 `Beta+/BA1+/XBB+`, `Beta+/BA+`) %>%
          rename(`BA.1 Only` = "BA1+", `XBB+` = "BA1+/XBB+",
                 `Proto+/Beta+/XBB+` = "Proto+/Beta+/BA1+/XBB+", `Proto+/Beta+` = "Proto+/Beta+/BA1+", `Proto+/XBB+` = "Proto+/BA1+/XBB+", `Proto+` = "Proto+/BA1+",
                 `Beta+/XBB+` = "Beta+/BA1+/XBB+", `Beta+` = "Beta+/BA+"
                 ) %>%
          mutate(TotalBA1 = rowSums(across(`BA.1 Only`:`Beta+`)),
                 across(c(`BA.1 Only`:`Beta+`), .fns = ~ . / TotalBA1)) %>%
  pivot_longer(!c(`Subject ID`, `Time point Guess`, Immunogen), names_to = "Population", values_to = "Proportion") %>%
  mutate(Population = factor(Population, levels = c(
    "BA.1 Only", "Proto+/Beta+/XBB+", "Proto+/Beta+", "Proto+/XBB+", "Proto+",
    "Beta+/XBB+", "Beta+", "XBB+"
  ))) %>%
  filter(Population != "TotalBA1") %>%
  arrange(`Time point Guess`) %>%
  group_by(Immunogen, `Subject ID`, Population) %>%
  mutate(FoldChange = Proportion / Proportion[[1]]) %>%
  filter(!is.nan(FoldChange), !is.infinite(FoldChange), !(`Time point Guess` %in% c(1, 90, 180))) %>%
  group_by(Immunogen, Population)%>%
  summarize(median = median(FoldChange))


ggplot(stats, aes(x= Population, y= Immunogen))+
  geom_point(aes(size = abs(median - 1), fill = median - 1), shape = 21)+
  scale_fill_gradient2(midpoint = 0, high = "#F05039", low = "#1F449C")+
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
ggsave(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure 2", "Figure2_BA1Binding_FoldChange.png"),width = 5, height = 2, units = "in", device = "png", dpi = 600)
dev.off()

#evaluate specific pops
###evaluate specific pops
stats <- flow %>% filter(probeset == "New") %>%
  select(`Subject ID`, Immunogen, `Time point Guess`, `BA1+`, `BA1+/XBB+`, 
         `Proto+/Beta+/BA1+/XBB+`, `Proto+/Beta+/BA1+`, `Proto+/BA1+/XBB+`, `Proto+/BA1+`,
         `Beta+/BA1+/XBB+`, `Beta+/BA+`) %>%
  rename(`BA.1 Only` = "BA1+", `XBB+` = "BA1+/XBB+",
         `Proto+/Beta+/XBB+` = "Proto+/Beta+/BA1+/XBB+", `Proto+/Beta+` = "Proto+/Beta+/BA1+", `Proto+/XBB+` = "Proto+/BA1+/XBB+", `Proto+` = "Proto+/BA1+",
         `Beta+/XBB+` = "Beta+/BA1+/XBB+", `Beta+` = "Beta+/BA+"
  ) %>%
  mutate(TotalBA1 = rowSums(across(`BA.1 Only`:`Beta+`)),
         across(c(`BA.1 Only`:`Beta+`), .fns = ~ . / TotalBA1)) %>%
  pivot_longer(!c(`Subject ID`, `Time point Guess`, Immunogen), names_to = "Population", values_to = "Proportion") %>%
  mutate(Population = factor(Population, levels = c(
    "BA.1 Only", "Proto+/Beta+/XBB+", "Proto+/Beta+", "Proto+/XBB+", "Proto+",
    "Beta+/XBB+", "Beta+", "XBB+"
  ))) %>%
  filter(Population != "TotalBA1") %>%
  arrange(`Time point Guess`) %>%
  group_by(Immunogen, `Subject ID`, Population) %>%
  mutate(FoldChange = Proportion / Proportion[[1]]) %>%
  filter(!is.nan(FoldChange), !is.infinite(FoldChange), !(`Time point Guess` %in% c(1, 90, 180)))

ggplot(stats[stats$FoldChange < 10,], aes(x = Population, y= FoldChange))+
  geom_boxplot(aes(fill = Population))+
  geom_point(aes(fill = Population), shape = 21, position = position_jitterdodge())+
  geom_hline(yintercept = 1)+
  facet_grid(cols = vars(Immunogen))+
  theme_classic()+
  theme(legend.key.size = unit(0.2, 'cm'),
        plot.title = element_text(size=6),
        axis.title.y = element_text(size=7),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=7,angle = 45, hjust=1.1, vjust=1.1),
        axis.text.y = element_text(size=7),
        strip.background = element_blank(),
        strip.text = element_text(size = 6, face = "bold", margin = margin()),
        panel.spacing = unit(0.3, "lines"),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 6),
        legend.box.spacing = margin(0.5))
