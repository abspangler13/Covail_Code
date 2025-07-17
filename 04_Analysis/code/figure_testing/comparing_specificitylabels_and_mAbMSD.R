library(ggplot2)
library(dplyr)
library(here)
library(Seurat)
library(readxl)
library(tidyseurat)
library(stringr)
library(alakazam)

set.seed(1)

#load in the data
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds"))

seuObj@meta.data <- seuObj@meta.data %>%
  mutate(PostExposure = case_when(Timepoint == "Day 15" ~ "Post-Boost",
                                  is.na(InfectionRange) ~ "Ignore",
                                  InfectionRange == "Between Days 15-90" & Timepoint == "Day 90" ~ "Post-Infection",
                                  InfectionRange == "Between Days 90-180" & Timepoint == "Day 180" ~ "Post-Infection",
                                  TRUE ~ "Ignore"))

#bring in the mAb MSD data- we specifically want to see if we can optimize the prototype-specific population
mAbData <- read.csv(here::here("04_Analysis", "data_objects", "figure_testing", "mAbData", "mAb_CITESeq_QC_Criteria.csv")) %>%
                filter(CiteSeq_Specificity == "Proto+Omi-")


#bring in the RATPIg data
#we want to also see if we can optimize this 
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

#####prepare sequencing+binding dataset
smartseqData <- read.csv(here::here("03a_Immcantation_SmartSeqMerge", "MergedCiteSeqAndSmartSeqSequences.csv")) %>%
  mutate(sample = paste0(str_extract(subject_id ,"[0-9]+"), "_", str_remove(Well_ID, "(?<=[A-H])0(?=[1-9])")),
         WellBarcode = sample)

#merge binding data into citeseq
short <- merged %>% select(c("BindingPopulation", "NeutralizingPopulations", "WellBarcode", "Prototype", "Omicron BA.1", "Omicron XBB", "JN1"))

overlap <- smartseqData %>%
  mutate(BindingPopulation = short$BindingPopulation[match(smartseqData$WellBarcode, short$WellBarcode)],
         ProbeBinding = case_when(!is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\+B\\+)|(Fully Cross-Reactive)") ~ "Proto+Omi+",
                                  !is.na(BindingPopulation) & str_detect(BindingPopulation, "(P\\-B\\+)|(BA.1-specific)") ~ "Proto-Omi+",
                                  is.na(BindingPopulation) & !is.na(adj.ProtoOmi) ~ adj.ProtoOmi,
                                  TRUE ~ "Smartseq Only")) %>%
  select(clone_subgroup_id, CELL, Dataset, ProbeBinding) %>%
  filter(ProbeBinding != "Smartseq Only") %>%
  group_by(clone_subgroup_id) %>%
  mutate(Specificity = case_when("SmartSeq" %in% Dataset ~ ProbeBinding[Dataset == "SmartSeq"][1],
                                 TRUE ~ ProbeBinding),
         InitialMistmatch = Specificity == ProbeBinding)

######
#time to tinker with the BA.1 cutoff and reattach the edited labels to cells in the smartseq data
#visualize populations
df <- seuObj@meta.data

ggplot(df, aes(x = `Proto-RBD-PE`+1, y = `BA1-RBD-PE`+1, fill = adj.ProtoOmi))+
  geom_point(size = 1, shape =21)+
  geom_hline(yintercept = 5+1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_classic()

df <- df %>% mutate(BAEdit = case_when(`BA1-RBD-PE` >= 4 ~ TRUE, TRUE ~ FALSE)) %>%
             mutate(edit.ProtoOmi = case_when((BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & `Proto.RBD.PE_Positive` ~ "Proto+Omi+",
                                             (!(BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & `Proto.RBD.PE_Positive` ~ "Proto+Omi-"),
                                             (BAEdit == TRUE | XBB.RBD.no.fluor_Positive == TRUE) & !`Proto.RBD.PE_Positive` ~ "Proto-Omi+",
                                             TRUE ~ "Proto-Omi-"))

adj.spec.P <- as.data.frame(df %>%
                              group_by(clone_subject_id,edit.ProtoOmi) %>% 
                              dplyr::summarise(Freq = n()) %>% 
                              pivot_wider(names_from = edit.ProtoOmi, values_from = Freq) %>% 
                              replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$clone_subject_id #this block takes the summary table and chooses the max label as the clonal specificity
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P$edited.ProtoOmi<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("edited.ProtoOmi","clone_subject_id")]
df2 <- df %>% left_join(adj.spec.P,by="clone_subject_id")

ggplot(df2, aes(x = `Proto-RBD-PE`+1, y = `BA1-RBD-PE`+1, fill = adj.ProtoOmi))+
  geom_point(size = 1, shape =21)+
  geom_hline(yintercept = 4+1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_classic()

ggplot(df2, aes(x = `Proto-RBD-PE`+1, y = `BA1-RBD-PE`+1, fill = edited.ProtoOmi))+
  geom_point(size = 1, shape =21)+
  geom_hline(yintercept = 4+1)+
  scale_y_log10()+
  scale_x_log10()+
  theme_classic()

ggplot(df2, aes(fill=adj.ProtoOmi, y=1, x=Booster, color=adj.ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of Probe Binding")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster- Old Labels")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))

ggplot(df2, aes(fill=edited.ProtoOmi, y=1, x=Booster, color=edited.ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Proportion of Probe Binding")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster- New Labels")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
####

#####
#does the new adj.Population label match the expected mAb binding profile?
mAbData$NewCiteseqLabel <- df2$edited.ProtoOmi[match(mAbData$cell, df2$CELL)]

mAbData2 <- mAbData %>% pivot_longer(!c(X, name, cell, V.D.J.REGION, Booster, Timepoint, c_call, CiteSeq_Specificity, ClusterLabel, ReasonForChoosing, NewCiteseqLabel),names_to = "Antigen", values_to = "ECL")

ggplot(mAbData2, aes(x = name, y = Antigen, fill = ECL))+
  geom_tile()+
  scale_y_discrete(limits = c("WA.1", "BA.1", "XBB"))+
  ggtitle("New Labels")+
  facet_grid(cols = vars(NewCiteseqLabel), scales = "free_x")+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
#####

#how does it stack up with flow?
#load in the flow data and make the necessary variables
flowRaw <- read_xlsx(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Unfiltered_COVAILFlowDataset_250318.xlsx")) %>% filter(Dataset != "Delta Panel")

flow <- flowRaw %>%
  mutate(TotalRBD = rowSums(select(.,contains("Combined"))),
         ProtoNotBeta = rowSums(select(., contains("Proto"), -contains("Beta"))),
         BetaNotProto = rowSums(select(., contains("Beta"), -contains("Proto"))),
         ProtoBeta = rowSums( select(.,matches("Proto.+Beta"))),
         ProtoNotOmicron = rowSums(select(., contains("Proto"), -contains("BA1"))),
         OmiNotProto = rowSums(select(., contains("BA1"), -contains("Proto"))),
         ProtoOmi = rowSums(select(.,matches("Proto.+BA"))),
         NotProto = rowSums(select(., c("Combined_Beta+/BA1+", "Combined_Beta+", "Combined_BA1+"))),
         PropProtoOmi = ProtoOmi / TotalRBD,
         Subject_Time = paste0(`Subject ID`, "_Day ", Timepoint)
         )

check <- df2 %>%
          group_by(Subject, Timepoint, edited.ProtoOmi) %>%
          summarize(n = n()) %>%
          mutate(Proportion = n / sum(n),
                 Subject_Time = case_when(Timepoint != "Day 0" ~ paste0(Subject, "_", Timepoint),
                                          TRUE ~ paste0(Subject, "_Day 1")))

flow2 <- flow %>%
          mutate(CITESeqProp = check$Proportion[match(.$Subject_Time, check$Subject_Time)]) %>%
  filter(!is.na(CITESeqProp))

ggplot(flow2, aes(x = PropProtoOmi, y = CITESeqProp, fill = Booster))+
  geom_point(shape =21)+
  geom_abline(slope = 1)+
  ylab("Proportion Cross-Reactive by CITESeq")+
  xlab("Proportion Cross-Reactive by Flow")+
  ylim(0,1)+
  xlim(0,1)+
  theme_classic()+
  theme()
