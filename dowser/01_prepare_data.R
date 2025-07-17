# https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html
# load libraries
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(sessioninfo))

seurat <- readRDS(file = here::here("dowser","covObj_clustered_demultiplexed.rds"))
seurat$cell_id <- Cells(seurat)

# df <- seurat@meta.data

# clones <- df[df$Infection == "N",] %>%
#             group_by(Booster, clone_subject_id, Timepoint) %>%
#             summarize(n = n()) %>%
#             mutate(Total = sum(n),
#                    TimepointCorrect = case_when("Day 180" %in% unique(Timepoint) & length(intersect(unique(Timepoint), c("Day 0", "Day 15", "Day 90"))) >= 1 ~ "Present",
#                                                 TRUE ~ "Nope"),
#                    Select = Total >= 15 & TimepointCorrect == "Present")

# check <- clones[clones$Select == TRUE,] #check outcome
# unique(check$clone_subject_id[check$Select == TRUE])
# # [1] "4957484948_993_207"  "5048574848_1658_61"  "5048574848_583_4"   
# # [4] "5048574848_976_188"  "4953494948_1968_130" "4953494948_2324_276"
# # [7] "4953494948_2720_226" "4953494948_3124_110" "4953494948_3191_51" 
# # [10] "4953494948_537_112"  "5357484948_2599_193" "5553564848_1584_37" 
# # [13] "4955534848_1795_127" "4955534848_514_95"   "4955534848_803_162" 
# # [16] "4955534848_850_19"

# my.data <- df[df$clone_subject_id %in% check$clone_subject_id,]

# my.clones = formatClones(my.data, clone = "clone_subject_id",traits ="Timepoint.num",text_fields = c("Subject"))

## load all heavy chain data# List all files matching the pattern
file_list <- list.files(path = here::here("..","..","Rory","Covail","03_Immcantation"), 
                                          pattern = "*heavy_light_clone-pass_germ-pass.tsv", 
                                          full.names = TRUE)

# Read all files and combine them into one dataframe
combined_df <- vroom::vroom(file_list, id = "filename") %>%
                                          bind_rows()
combined_df$filename <- basename(combined_df$filename)
combined_df$subject_id <- sub("_heavy_light_clone-pass_germ-pass.tsv", "", combined_df$filename)

saveRDS(combined_df, file = here::here("dowser","combined_heavy_light_clone-pass_germ-pass.rds"))

## load all light chain data# List all files matching the pattern
file_list <- list.files(path = here::here("..","..","Rory","Covail","03_Immcantation"), 
                                          pattern = "*light_parse-select.tsv", 
                                          full.names = TRUE)

# Read all files and combine them into one dataframe
lc_combined_df <- vroom::vroom(file_list, id = "filename") %>%
                                          bind_rows()
lc_combined_df$filename <- basename(lc_combined_df$filename)
lc_combined_df$subject_id <- sub("_light_parse-select.tsv", "", lc_combined_df$filename)

saveRDS(lc_combined_df, file = here::here("dowser","combined_light_parse-select.rds"))

# Combine the heavy and light chain dataframes
combined_all_df <- bind_rows(combined_df, lc_combined_df)

# Save the combined dataframe
saveRDS(combined_all_df, file = here::here("dowser","all_bcr_data.rds"))
# combined_all_df <- readRDS(file = here::here("dowser","all_bcr_data.rds"))

# Filter bcr data to only contain the sequences that come from cells that are contained in the seruat object
filtered_combined_all_df <- combined_all_df %>% 
    filter(cell_id %in% seurat$cell_id)

dim(filtered_combined_all_df)
#[1] 40320    43

### Add meta data
# Add these columns from seurat@meta.data to filtered_combined_all_df based on cell_id
meta_data_columns <- c("Population", "adj.Population","ProtoOmi", "adj.ProtoOmi", "TimeAdj.ProtoOmi", "harmony.snn_res.0.4", 
                       "seurat_clusters", "SecondClusteringLabel", "InfectionTimepoint", "InfectionRange", 
                       "ClusterLabel", "EnrollmentDate", "InfectionDate", "InfectionLineage", "InfectionTrunc", 
                       "InfectionIndentMethod", "MULTI_classification", "Timepoint", "Subject", "Booster", 
                       "Infection", "Bcell","TimeAdj.ProtoOmi")

filtered_combined_all_df <- filtered_combined_all_df %>%
    left_join(seurat@meta.data[, c("cell_id", meta_data_columns)], by = "cell_id")

filtered_combined_all_df <-filtered_combined_all_df %>% mutate(Timepoint.num = case_when(Timepoint == "Day 0" ~ 0,
                                                        Timepoint == "Day 15" ~ 15,
                                                        Timepoint == "Day 90" ~ 90,
                                                        Timepoint == "Day 180" ~ 180))

# change clone_id to contain subject id
filtered_combined_all_df$clone_id_orig <- filtered_combined_all_df$clone_id
filtered_combined_all_df$clone_id <- paste0(filtered_combined_all_df$subject_id,"_",filtered_combined_all_df$clone_id)

saveRDS(filtered_combined_all_df, file = here::here("dowser","all_bcr_data_meta.rds"))
# filtered_combined_all_df <- readRDS(file = here::here("dowser","all_bcr_data_meta.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()