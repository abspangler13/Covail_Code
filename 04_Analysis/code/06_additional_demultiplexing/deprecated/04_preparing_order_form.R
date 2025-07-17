#prepare the sequences
library(here)
library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(Seurat)

#load in the imgt data
hc <- read_xls(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "imgt", "vquest_hc.xls"))
lc <- read_xls(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "imgt", "vquest_lc.xls"))

#load in seurat object so we can get info on stuff (i should've done this before woopsies)
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
df <- seuObj@meta.data

#create a column that specifies cells for hc and lc (remove the _HC and _IGK/IGL tags)
hc$cell <- str_extract(hc$`Sequence ID`, ".*(?= _HC|_IGL|_IGK)") #should all be HC but added LC tags for ease of mind
lc$cell <- str_extract(lc$`Sequence ID`, ".*(?= _HC|_IGL|_IGK)")

#append true names to each chain- "CoV-Donor-Clone"
hc$name <- paste0("CoV_", df$Subject[match(hc$cell, df$CELL)], "_", df$clone_id[match(hc$cell, df$CELL)], "_IGH")
lc$name <- paste0("CoV_", df$Subject[match(lc$cell, df$CELL)], "_", df$clone_id[match(lc$cell, df$CELL)], 
                  ifelse(str_detect(lc$`Sequence ID`,"IGL"), "_IGL", "_IGK"))
lc$CorrectedVDJ <- ifelse(str_detect(lc$`Sequence ID`,"IGK"), lc$`V-J-REGION`, paste0(lc$`V-J-REGION`, "GQPKANPTVTLFP"))

#add to output dataframe
lcSmall <- lc %>% select(name, cell, CorrectedVDJ) %>% rename(`V-D-J-REGION` = CorrectedVDJ)

output <- hc %>% select(name, cell, `V-D-J-REGION`) %>% rbind(lcSmall)

#save
write_xlsx(output, here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "COVAIL_mAbOrdering_Information.xlsx"))
