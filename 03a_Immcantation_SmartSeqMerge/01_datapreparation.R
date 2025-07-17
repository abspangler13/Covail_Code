#load dependencies
library(dplyr)
library(readxl)
library(here)
library(stringr)
library(scoper)
library(alakazam)
library(dowser)

########load in the smartseq data
#smartseq dataset 1
#smartseq1 <- read.csv(here::here("01_raw-data", "RATPIg", "SmartSeq", "simplified_bcr_data 2.csv"))
smartseq1 <- read.csv(here::here("01_raw-data", "RATPIg", "SmartSeq", "bcr_data_top_heavy_light_2.csv"))

#smartseq dataset 2
#smartseq2 <- read_xlsx(here::here("01_raw-data", "RATPIg", "SmartSeq", "simplified_bcr_data 3.xlsx"))
smartseq2 <- read.csv(here::here("01_raw-data", "RATPIg", "SmartSeq", "bcr_data_top_heavy_light_run_3.csv"))

#combine the two datasets
df <- rbind(smartseq1, smartseq2) %>%
        mutate(Well_ID = substr(MiXCR_tagValueCELL0ROW,start = 1,stop = 3),
               donor_id = case_when(subject_id == "COVAIL_205772412" ~ "4951574848",
                                    subject_id == "COVAIL_206230717" ~ "5752574848",
                                    subject_id == "COVAIL_206231302" ~ "5657574848",
                                    subject_id == "COVAIL_206268325" ~ "5557544848",
                                    subject_id == "COVAIL_206311727" ~ "5249544848",
                                    subject_id == "COVAIL_206319448" ~ "4955534848",
                                    subject_id == "COVAIL_206333104" ~ "5356534848",
                                    subject_id == "COVAIL_206355391" ~ "5048574848",
                                    subject_id == "COVAIL_206355448" ~ "4848544848",
                                    subject_id == "COVAIL_206362879" ~ "5556484948",
                                    subject_id == "COVAIL_206363191" ~ "5553554848",
                                    subject_id == "COVAIL_206363266" ~ "4957484948",
                                    subject_id == "COVAIL_206372353" ~ "4957564848",
                                    subject_id == "COVAIL_206382618" ~ "5755544848",
                                    subject_id == "COVAIL_206383247" ~ "5750564848",
                                    subject_id == "COVAIL_206383617" ~ "5050484948",
                                    subject_id == "COVAIL_206383965" ~ "5457484948",
                                    subject_id == "COVAIL_206447647" ~ "5048544848",
                                    subject_id == "COVAIL_205600941top_206292415bottom" & str_detect(Well_ID, "[A-D]") ~ "4954554848",
                                    subject_id == "COVAIL_205600941top_206292415bottom" & str_detect(Well_ID, "[E-H]") ~ "5351564848",
                                    TRUE ~ "Woops"),
               CELL = paste(donor_id, Well_ID, sep="_"))

#correct the well ID for the 206362879 plate! It was incorrectly transfected (woops)
errorDonor <- df %>% filter(subject_id == "COVAIL_206362879") %>%
                mutate(row  = str_extract(Well_ID, "[A-H]"),
                       column = as.numeric(str_extract(Well_ID, "[0-9]+")),
                       row = case_when(column == 12 & row == "A" ~ "H",
                                       column == 12 & row == "B" ~ "G",
                                       column == 12 & row == "C" ~ "F",
                                       column == 12 & row == "D" ~ "E",
                                       column == 12 & row == "E" ~ "D",
                                       column == 12 & row == "F" ~ "C",
                                       column == 12 & row == "G" ~ "B",
                                       column == 12 & row == "H" ~ "A",
                                       TRUE ~ row),
                       column = case_when(column %in% c(1,2,3,4,5,6,7,8,9,10) ~ column + 1,
                                          column == 12 ~ 12,
                                          column == 11 ~ 1),
                       column = as.character(column),
                       column = case_when(column %in% c("0","1","2","3","4","5","6","7","8","9") ~ paste0("0", column),
                                          TRUE ~ column),
                       Well_ID = paste0(row, column)) %>%
                    select(!c("row", "column"))

#merge updated well_id from error plate
mergeddf <- df %>% filter(subject_id != "COVAIL_206362879") %>%
  rbind(errorDonor)

########load in the CITESeq data- I'm using this instead of the raw tsvs as we only want clonality calls for cells we know are in our data
seuObj <- readRDS(file=here::here("04_Analysis", "data_objects", "04_probe", "CoVSeuratObj_VDJCSOGEX_SpecificitiesLabelled_CloneCorrected.rds"))
CITESeqData <- seuObj@meta.data %>% mutate(donor_id = Subject) %>% select(intersect(c(colnames(df), paste0("LC_",colnames(df))), colnames(seuObj@meta.data)), donor_id, CELL, adj.ProtoOmi)

# mergedDF <- bind_rows(df, CITESeqData)

########we need to pivot light chain and heavy chain information to a longer format
light_chain <- CITESeqData %>% select("CELL", "donor_id", contains("LC_"))
heavy_chain <- CITESeqData %>% select("CELL", "donor_id", !contains("LC_"))

colnames(light_chain) <- str_remove(colnames(light_chain), "LC_")

#bind rows
mergedPivoted <- bind_rows(heavy_chain, light_chain, mergeddf)

########perform clonality analysis using scoper
results <- hierarchicalClones(mergedPivoted,
                              cell_id = "CELL",
                              method = "nt",
                              threshold = 0.15,
                              only_heavy = TRUE, split_light = TRUE,
                              summarize_clones = FALSE,
                              fields = "donor_id")

#split by light chain- adds sub clone id
results$v_gene_call = alakazam::getGene(results$v_call)
results$j_gene_call = alakazam::getGene(results$j_call)
split = resolveLightChains(results, j_call = "j_gene_call", v_call = "v_gene_call", cell = "CELL")


########okay slay it worked!! recombine the data in a wide format
light_chain_pivot <- split %>% filter(str_detect(v_call, "K|L"))

wideData <- split %>% filter(!str_detect(v_call, "K|L")) %>%
            merge(., light_chain_pivot, by = "CELL", suffixes = c("","_LC"))

########label dataset origin
wideData$Dataset <- ifelse(str_detect(wideData$sequence_id, "COVAIL"), "SmartSeq", "CITESeq")

########write the data
write.csv(wideData, here::here("03a_Immcantation_SmartSeqMerge", "MergedCiteSeqAndSmartSeqSequences.csv"))
