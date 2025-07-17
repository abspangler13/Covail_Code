library(dplyr)
library(here)
library(tidyverse)
library(readxl)
library(writexl)
library(stringr)

myRawData <- read.csv(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Unfiltered_COVAILFlowDataset_250318.csv"))
myFilteredData <- read.csv(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "Filtered_COVAILDataset_250318_CORRECTED.csv"))

sarahData <- read_xlsx(here::here("01_raw-data", "FlowData", "FinalizedDatasets", "sarahData.xlsx"))

#people in mydata but not sarahdata
myRawDataSpecific <- setdiff(myRawData$Specimen.ID, sarahData$`Specimen ID`)
myFilteredDataSpecific <- setdiff(myFilteredData$Specimen.ID, sarahData$`Specimen ID`)
filterVsNot <- setdiff(myRawData$Specimen.ID, myFilteredData$Specimen.ID)
sarahDataSpecificRaw <- setdiff(sarahData$`Specimen ID`, myRawData$Specimen.ID)
sarahDataSpecificFiltered <- setdiff(sarahData$`Specimen ID`, myFilteredData$Specimen.ID)

#add necessary variables used in filtering step
myRawData$RepeatBarcode <- paste(myRawData$Subject.ID, myRawData$Timepoint)
repeats <- myRawData$Subject.ID[duplicated(myRawData$RepeatBarcode)]

#Label each step that data would be filtered out
oosBoost <- unique(myRawData$Subject.ID[myRawData$oosboost_flag == "Y" & !is.na(myRawData$infect_flag)])
infect <- unique(myRawData$Subject.ID[myRawData$infect_flag == "Y" & !is.na(myRawData$infect_flag)])
day0 <- myRawData$Subject.ID[myRawData$Timepoint == 1]
day15 <- myRawData$Subject.ID[myRawData$Timepoint == 15]
missing <- myRawData$Subject.ID[!(myRawData$Subject.ID %in% day0 & myRawData$Subject.ID %in% day15)]

myRawData <- myRawData %>%
              mutate(Outcome = case_when(
                is.na(Specimen.ID) ~ "Specimen ID is NA",
                Subject.ID %in% oosBoost ~ "Subject received out of study boost",
                Subject.ID == 5752524948 ~ "Issue with memory titers- very few MBCs",
                Subject.ID == 5349564848 & Timepoint == "15" ~ "Issue with viability",
                Subject.ID %in% repeats & Dataset == "Delta Panel" ~ "Repeat b/w XBB and Delta; XBB kept",
                infect_baseline == "Y" ~ "Infection at baseline",
                infect_flag == "Y" ~ "Infect_flag is 'Y'",
                Subject.ID %in% missing ~ "Subject does not have both day 0 and day 15 data",
                is.na(Timepoint) ~ "Missing timepoint data",
                Timepoint == "71" ~ "Day 71 timepoint",
                Booster == "Delta + BA.1 mRNA" ~ "Delta + BA.1 mRNA",
                Subject.ID %in% c(5353564848, 5450574848) ~ "Donor has day 1 data for Delta panel but not XBB panel; run on both for only some timepoints",
                TRUE ~ "Included in filtered data"
              ),
              SarahDataset = case_when(
                Specimen.ID %in% myRawDataSpecific ~ "In my raw data but not Sarah's",
                Specimen.ID %in% filterVsNot ~ "Present in raw data and Sarah's, but not the filtered data",
                TRUE ~ "Present in all datasets"
              ))

write.csv(myRawData, here::here("01_raw-data", "FlowData", "ComparingSarahsDatasetAndRorys.csv"))

write.csv(table(myRawData$Outcome, myRawData$SarahDataset), here::here("01_raw-data", "FlowData", "ComparingSarahsDatasetAndRorys_table.csv"))

#sarah wants specimen IDs for rows 6 and 7
check <- myRawData %>% filter(Outcome %in% c("Infect_flag is 'Y'", "Subject does not have both day 0 and day 15 data")) %>%
          select(Specimen.ID, Subject.ID, Dataset, Outcome)

write.csv(check, here::here("01_raw-data", "FlowData", "QC", "SpecimenIDs_rows6and7.csv"))


#check
sarahData$SarahDataset <- case_when(sarahData$`Specimen ID` %in% sarahDataSpecificRaw ~ "Not in my data",
                                    TRUE ~ "In my data")


#check sarah's data- she says she does actually have timepoint info?
check2 <- sarahData %>% filter(`Specimen ID` %in% myRawData$Specimen.ID[myRawData$Outcome == "Missing timepoint data"])
