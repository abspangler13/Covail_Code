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

#####prepare sequencing+binding dataset
smartseqData <- read.csv(here::here("04_Analysis", "data_objects", "figure_testing", "testing_ratpig_inclusion", "MergedCiteSeqAndSmartSeqSequences_rbdneg.csv")) %>%
  mutate(sample = paste0(str_extract(subject_id ,"[0-9]+"), "_", str_remove(Well_ID, "(?<=[A-H])0(?=[1-9])")),
         WellBarcode = sample)

#merge binding data into citeseq
combinedData <- smartseqData %>% left_join(., merged, by = "WellBarcode")

#add in IgG data
igg <- read_xlsx(here::here("01_raw-data", "RATPIg", "COVAIL_Interpolated_IgGMSD_data.xlsx")) %>% filter(!is.na(Dilution_fold)) %>%
          mutate(WellBarcode = paste0(Donor, "_", WellID))

#merge into the combined data
overlapIgg <- combinedData %>%
                mutate(IgGQuant = case_when(Dataset == "CITESeq" ~ NA,
                                            TRUE ~ igg$CalculatedConcentration[match(.$WellBarcode, igg$WellBarcode)]))

#write xlsx of CITESeq + MSD + neuts
write_xlsx(overlapIgg, here::here("04_Analysis", "data_objects", "figure_testing", "testing_ratpig_inclusion", "CombinedMSD_Neut_Smartseq_OverlappingCITESeq.xlsx"))
#####