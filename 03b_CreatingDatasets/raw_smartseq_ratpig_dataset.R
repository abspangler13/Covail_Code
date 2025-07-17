library(ggplot2)
library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

#####read in MSD data
df <- readxl::read_xlsx(here::here("01_raw-data", "RATPIg", "MSD", "250428_CombinedMSDdata_KJG.xlsx")) |>
  filter(!str_detect(`RATPIg Well`, "(blank)|(VRC01)|(Positive)")) %>%
  mutate(CorrectedBoost = case_when(`Boost Type` %in% c("Proto+BA1", "Prototype+BA.1") ~ "BA.1",
                                    TRUE ~ `Boost Type`),
         `Omicron BA.1` = ifelse(`Omicron BA.1.old` < `Omicron BA.1.new`, `Omicron BA.1.old`, `Omicron BA.1.new`))

threshold = 20000

filtered_df <- df |>
  filter(!(`Donor ID` == 205772412 & `MSD Date` == as.Date("2024-12-20"))) |> #remove a donor who we ran twice
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

#add in IgG data
igg <- read_xlsx(here::here("01_raw-data", "RATPIg", "COVAIL_Interpolated_IgGMSD_data.xlsx")) %>% filter(!is.na(Dilution_fold)) %>%
  mutate(WellBarcode = paste0(Donor, "_", WellID))

#add in igg values
merged$IgGQuant <- igg$CalculatedConcentration[match(merged$WellBarcode, igg$WellBarcode)]

####
#####load in sequencing data
smartseqData <- read.csv(here::here("03a_Immcantation_SmartSeqMerge", "MergedCiteSeqAndSmartSeqSequences.csv")) %>% filter(Dataset == "SmartSeq") %>%
  mutate(subject_id = case_when(
    str_detect(CELL, "5351564848_E") ~ str_replace(CELL, "5351564848_E", "206292415_A"),
    str_detect(CELL, "5351564848_F") ~ str_replace(CELL, "5351564848_F", "206292415_B"),
    str_detect(CELL, "5351564848_G") ~ str_replace(CELL, "5351564848_G", "206292415_C"),
    str_detect(CELL, "5351564848_H") ~ str_replace(CELL, "5351564848_H", "206292415_D"),
    TRUE ~ subject_id
  ),
  Well_ID = case_when(
    str_detect(CELL, "5351564848_E") ~ str_replace(Well_ID, "E", "A"),
    str_detect(CELL, "5351564848_F") ~ str_replace(Well_ID, "F", "B"),
    str_detect(CELL, "5351564848_G") ~ str_replace(Well_ID, "G", "C"),
    str_detect(CELL, "5351564848_H") ~ str_replace(Well_ID, "H", "D"),
    TRUE ~ Well_ID
  ),
    #we need to fix the partial plate labels that I combined
        sample = paste0(str_extract(subject_id ,"[0-9]+"), "_", str_remove(Well_ID, "(?<=[A-H])0(?=[1-9])")),
         WellBarcode = sample)

#create a column describing presence in SmartSeq data
updated <- merged %>% 
            mutate(SmartseqPresent = WellBarcode %in% smartseqData$WellBarcode,
                   RemovedFromFilteredData = `Omicron BA.1` <= threshold)

#write data
write_xlsx(updated, here::here("01_raw-data", "COVAIL_unfiltered_MSD_Neut_IgG_Data.xlsx"))
#####