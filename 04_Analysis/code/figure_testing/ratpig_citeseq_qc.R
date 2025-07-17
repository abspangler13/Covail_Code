library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

#####read in MSD data
df <- readxl::read_xlsx(here::here("01_raw-data", "RATPIg", "MSD", "COVAIL_RATPIgMSD_FullDataSet_KJG_RM.xlsx"), sheet = "Data") |>
  filter(!str_detect(`RATPIg Well`, "(blank)|(VRC01)|(Positive)")) %>%
  mutate(CorrectedBoost = case_when(`Boost Type` %in% c("Proto+BA1", "Prototype+BA.1") ~ "BA.1",
                                    TRUE ~ `Boost Type`))

threshold = 20000

filtered_df <- df |>
  filter(`Omicron BA.1`>= threshold, !(`Donor ID` == 205772412 & `MSD Date` == as.Date("2024-12-20"))) |>
  mutate(`Donor ID` = as.character(`Donor ID`))

labels <- filtered_df |>
  mutate(BindingPopulation = case_when((`Omicron BA.1` >= threshold) & (Prototype >= `Omicron BA.1`*0.2) & (`Omicron XBB` >= `Omicron BA.1`*0.2) & (`JN1` >= `Omicron BA.1`*0.2) ~ "Fully Cross-Reactive",
                                       `Omicron BA.1` >= threshold & Prototype >= `Omicron BA.1`*0.2 & `Omicron XBB` >= `Omicron BA.1`*0.2 ~ "P+B+X+J-",
                                       `Omicron BA.1` >= threshold & Prototype >= `Omicron BA.1`*0.2 & `JN1` >= `Omicron BA.1`*0.2 ~ "P+B+X-J+",
                                       `Omicron BA.1` >= threshold & Prototype >= `Omicron BA.1`*0.2 ~ "P+B+X-J-",
                                       `Omicron XBB` >= `Omicron BA.1`*0.2 & `JN1` >= 0.2*`Omicron BA.1` ~ "P-B+X+J+",
                                       `Omicron XBB` >= `Omicron BA.1`*0.2  ~ "P-B+X+J-",
                                       `JN1` >= 0.2*`Omicron BA.1` ~ "P-B+X-J+",
                                       `Omicron BA.1` >= threshold ~ "BA.1-specific",
                                       TRUE ~ "Woopsies"),
         BindingPopulation = factor(BindingPopulation, levels = c("Fully Cross-Reactive", "P+B+X+J-", "P+B+X-J+", "P+B+X-J-", "P-B+X+J+", "P-B+X+J-", "P-B+X-J+", "BA.1-specific")),
         WellBarcode = paste0(`Donor ID`,"_", `RATPIg Well`))

#read in and append the neutralization data
neutralization <- read_xlsx(here::here("01_raw-data", "RATPIg", "Neut", "COVAIL RATPIG SUP SARSNeut_Binding_IgG DataUpdated20250205.xlsx")) %>%
                    mutate(`SampleSourcePlateWell ID` = str_remove(`SampleSourcePlateWell ID`,"(?<=[A-H])0"),
                           WellBarcode = paste0(SampleSourcePlateID, "_", `SampleSourcePlateWell ID`))

#merge the MSD and neut data
thresholdN <- 0.5

merged <- left_join(labels, neutralization, by = "WellBarcode") %>%
            mutate(NeutralizingPopulations = case_when(D614G >= thresholdN & BA.1 >= thresholdN & `XBB.1.5` >= thresholdN & `JN.1` >= thresholdN ~ "P+B+X+J+",
                                                       D614G >= thresholdN & BA.1 >= thresholdN & `XBB.1.5` >= thresholdN ~ "P+B+X+J-",
                                                       D614G >= thresholdN & BA.1 >= thresholdN & `JN.1` >= thresholdN ~ "P+B+X-J+",
                                                       D614G >= thresholdN & BA.1 >= thresholdN ~ "P+B+X-J-",
                                                       D614G >= thresholdN ~ "P+B-X-J-",
                                                       BA.1 >= thresholdN & `XBB.1.5` >= thresholdN & `JN.1` >= thresholdN ~ "P-B+X+J+",
                                                       BA.1 >= thresholdN & `XBB.1.5` >= thresholdN ~ "P-B+X+J-",
                                                       BA.1 >= thresholdN & `JN.1` >= thresholdN ~ "P-B+X-J+",
                                                       BA.1 >= thresholdN ~ "P-B+X-J-",
                                                       `XBB.1.5` >= thresholdN & `JN.1` >= thresholdN ~ "P-B-X+J+",
                                                       `XBB.1.5` >= thresholdN ~ "P-B-X+J-",
                                                       `JN.1` >= thresholdN ~ "P-B-X-J+",
                                                       TRUE ~ "Non-Neutralizing"))

# #write dataset
# write.csv(merged, here::here("04_Analysis", "data_objects", "paperfigures", "Figure 5", "COVAIL_MergedBinding_And_Neutralization.csv"))

#read in citeseq data
citeseq <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))

df2 <- citeseq@meta.data

# ggplot(df, aes(x = `Proto-RBD-PE`+1, y = `BA1-RBD-PE`+1, fill = adj.ProtoOmi))+
#   geom_point(size = 1, shape =21)+
#   geom_hline(yintercept = 5+1)+
#   scale_y_log10()+
#   scale_x_log10()+
#   theme_classic()
# 
# df <- df %>% mutate(BAEdit = case_when(`BA1-RBD-PE` >= 4 ~ TRUE, TRUE ~ FALSE)) %>%
#   mutate(edit.ProtoOmi = case_when((BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & `Proto-RBD-PE` ~ "Proto+Omi+",
#                                    (!(BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & `Proto-RBD-PE` ~ "Proto+Omi-"),
#                                    (BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & !`Proto-RBD-PE` ~ "Proto-Omi+",
#                                    TRUE ~ "Proto-Omi-"))
# 
# adj.spec.P <- as.data.frame(df %>%
#                               group_by(clone_subject_id,edit.ProtoOmi) %>% 
#                               dplyr::summarise(Freq = n()) %>% 
#                               pivot_wider(names_from = edit.ProtoOmi, values_from = Freq) %>% 
#                               replace(is.na(.),0))
# 
# rownames(adj.spec.P) <- adj.spec.P$clone_subject_id #this block takes the summary table and chooses the max label as the clonal specificity
# adj.spec.P <- adj.spec.P[,-1]
# adj.spec.P$adj.ProtoOmi<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
# adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
# adj.spec.P <- adj.spec.P[,c("adj.ProtoOmi","clone_subject_id")]
# df2 <- df %>% left_join(adj.spec.P,by="clone_subject_id")
# #####

###################Workspace
######
#read in smartseq data (includes overlaps with CITESeq) and verify probe+ labels
#####
smartseqData <- read.csv(here::here("03a_Immcantation_SmartSeqMerge", "MergedCiteSeqAndSmartSeqSequences.csv")) %>%
                mutate(sample = paste0(str_extract(subject_id ,"[0-9]+"), "_", str_remove(Well_ID, "(?<=[A-H])0(?=[1-9])")),
                       WellBarcode = sample)

#merge binding data into citeseq
short <- merged %>% select(c("BindingPopulation", "NeutralizingPopulations", "WellBarcode", "Prototype", "Omicron BA.1", "Omicron XBB", "JN1"))

seqProt <- smartseqData %>%
            mutate(BindingPopulation = short$BindingPopulation[match(smartseqData$WellBarcode, short$WellBarcode)],
                   NeutralizingPopulations = short$NeutralizingPopulations[match(smartseqData$WellBarcode, short$WellBarcode)],
                   Proto = short$Prototype[match(smartseqData$WellBarcode, short$WellBarcode)],
                   BA1 = short$`Omicron BA.1`[match(smartseqData$WellBarcode, short$WellBarcode)],
                   XBB = short$`Omicron XBB`[match(smartseqData$WellBarcode, short$WellBarcode)],
                   adj.ProtoOmi = df2$adj.ProtoOmi[match(smartseqData$CELL, df2$CELL)])

#how many overlapping clonal groups get the probe-binding correct?
stats <- seqProt %>%
          select(Dataset, clone_subgroup_id, adj.ProtoOmi, BindingPopulation) %>%
          mutate(ProbeBinding = case_when(!is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\+B\\+)|(Fully Cross-Reactive)") ~ "Proto+Omi+",
                                          !is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\-B\\+)|(BA.1-specific)") ~ "Proto-Omi+",
                                          is.na(BindingPopulation) & !is.na(adj.ProtoOmi) ~ adj.ProtoOmi,
                                          TRUE ~ "Smartseq Only")) %>%
          filter(ProbeBinding != "Smartseq Only") %>%
          group_by(clone_subgroup_id, Dataset, ProbeBinding) %>%
          summarize(n = n()) %>%
          group_by(clone_subgroup_id) %>%
          mutate(Overlap = length(unique(Dataset)) > 1,
                 Disagreement = length(unique(ProbeBinding)) > 1) %>%
          filter(Overlap == TRUE)

s1 <- stats %>%
      group_by(Disagreement) %>%
      summarize(n = sum(n)) %>%
      ungroup() %>%
      mutate(prop = n / sum(n),
             Overlap = "Overlap",
             Disagreement = factor(Disagreement, levels = c("TRUE", "FALSE")))

#plot barplot of disagreement
ggplot(s1, aes(x = "Overlap", y = prop, fill = Disagreement))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values = c("FALSE" = "#9d02d7", "TRUE" = "#fa8775"))+
  ylab("Proportion")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position = "left")

#what are the most common types of disagreement?
s2 <- stats %>%
      filter(Disagreement == TRUE) %>%
      mutate(Identifier = paste0(Dataset,"_",ProbeBinding),
             Conflict = case_when(sum(Dataset == "CITESeq") > 1 ~ "CITESeq Labels Conflict",
                                  "CITESeq_Proto+Omi-" %in% Identifier & "SmartSeq_Proto+Omi+" %in% Identifier ~ "CITESeq BA.1-",
                                  "CITESeq_Proto+Omi+" %in% Identifier & "SmartSeq_Proto-Omi+" %in% Identifier ~ "CITESeq Prototype+",
             TRUE ~ "Error")) %>%
      group_by(clone_subgroup_id, Conflict) %>%
      summarize(n = n()) %>%
      group_by(Conflict) %>%
      summarize(n = n()) %>%
      mutate(Proportion = n / sum(n),
             Overlap = "Overlap")
      
ggplot(s2, aes(x = "Overlap", y = Proportion, fill = Conflict))+
  geom_bar(stat = "identity", position = "stack", color = "black")+
  scale_fill_manual(values = c("CITESeq BA.1-" = "#00429d",
                               "CITESeq Prototype+" = "#ffffe0",
                               "CITESeq Labels Conflict" = "red"))+
  ylab("Proportion")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        legend.position = "right",
        axis.text.x = element_blank())
####


#######
#load in citeseq data and make correction to BA.1 signal

