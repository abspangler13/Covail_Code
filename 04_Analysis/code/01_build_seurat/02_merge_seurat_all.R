library(Seurat)
library(here)
library(tidyverse)
library(gridExtra)
library(sessioninfo)    
library(tidyseurat)

#add run info
run_number <- 1
run_info <- data.frame(c("Run" = run_number))
# run_info <- read.csv(file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"))

#run names
sample_id = c("COV_09_13_2023_P1",
              "COV_09_19_2023_P1",
              "COV_09_19_2023_P2",
              "COV_09_19_2023_P3",
              "COV_09_20_2023_P1",
              "COV_09_20_2023_P2",
              "COV_09_20_2023_P3",
              "COV_09_21_2023_P1",
              "COV_09_21_2023_P2",
              "COV_09_26_2023_P1",
              "COV_09_26_2023_P2",
              "COV_09_27_2023_P1",
              "COV_09_27_2023_P2",
              "COV_09_28_2023_P1",
              "COV_09_28_2023_P2",
              "COV_09_28_2023_P3",
              "COV_10_05_2023_P1",
              "COV_10_06_2023_P1",
              "COV_10_06_2023_P2"
)

#for testing
#sample_id <- sample_id[1:2]

#Seurat object list
seuratPosObjs <- list()
seuratNegObjs <- list()

######## Read in positive RDS objects
for(i in 1:length(sample_id)){
  seuratPosObjs[i] <- readRDS(file = here::here("04_Analysis", "data_objects", "01_build_seurat", sample_id[i], paste0(sample_id[i], "_pos.rds")))
}

#set object names
for(i in 1:length(sample_id)){
  names(seuratPosObjs)[i] <- sample_id[i]
}

######## Read in negative RDS objects
for(i in 1:length(sample_id)){
  seuratNegObjs[i] <- readRDS(file = here::here("04_Analysis", "data_objects", "01_build_seurat", sample_id[i], paste0(sample_id[i], "_neg.rds")))
}

for(i in 1:length(sample_id)){
  names(seuratNegObjs)[i] <- sample_id[i]
}

#renames probes before merge- I've kept one chunk of pre-merge cleaning from Abby's code just to
#have a reference in case there are any issues- we used the exact same names/reagents for every run
#with the exception of some hashtags for a few different pools, so this could present an issue- these are noted on the previous script and excluded from the HTO csv
#P310.p
# HA.df <- GetAssayData(P310.n, assay = "Probes", slot = "counts")
# rownames(HA.df)
# rownames(HA.df)<-c("H2-oligo","SA-PE-0951","SA-0972","SA-PE-0953","SA-PE-0954","SA-PE-0955")
# rownames(HA.df)
# HA.assay<- CreateAssayObject(counts = HA.df)
# rownames(HA.assay)
# P310.n[["HAs"]] <- HA.assay
# P310.n[['Probes']] <- NULL
# P310.n <- RenameAssays(P310.n, HAs = "Probes")

#merge objects into one and save
#you can't coerce multiple items from a list using double brackets, which is the only way to
#correctly access the data. As a result, I'll need to input it myself
mergedPosObjs <- sp::merge(seuratPosObjs[[1]], y = c(seuratPosObjs[[2]], seuratPosObjs[[3]],
                                                     seuratPosObjs[[4]], seuratPosObjs[[5]],
                                                     seuratPosObjs[[6]], seuratPosObjs[[7]],
                                                     seuratPosObjs[[8]], seuratPosObjs[[9]],
                                                     seuratPosObjs[[10]], seuratPosObjs[[11]],
                                                     seuratPosObjs[[12]], seuratPosObjs[[13]],
                                                     seuratPosObjs[[14]], seuratPosObjs[[15]],
                                                     seuratPosObjs[[16]], seuratPosObjs[[17]],
                                                     seuratPosObjs[[18]], seuratPosObjs[[19]]), add.cell.ids = c(sample_id))

message("Counts for each pool in merged positive Seurat object: \n",paste0(capture.output(plyr::count(mergedPosObjs@meta.data$orig.ident)), collapse = "\n"))
#       x freq
# 1  M310 1408
# 2 M310B 2952
# 3  P310 1457
# 4 P310B  741

#merge into one object and save for negatives
mergedNegObjs <- merge(seuratNegObjs[[1]], y = c(seuratNegObjs[[2]], seuratNegObjs[[3]],
                                                 seuratNegObjs[[4]], seuratNegObjs[[5]],
                                                 seuratNegObjs[[6]], seuratNegObjs[[7]],
                                                 seuratNegObjs[[8]], seuratNegObjs[[9]],
                                                 seuratNegObjs[[10]], seuratNegObjs[[11]],
                                                 seuratNegObjs[[12]], seuratNegObjs[[13]],
                                                 seuratNegObjs[[14]], seuratNegObjs[[15]],
                                                 seuratNegObjs[[16]], seuratNegObjs[[17]],
                                                 seuratNegObjs[[18]], seuratNegObjs[[19]]), add.cell.ids = c(sample_id))

#we can add timepoints to our objects now! every hashtag defines the same timepoint regardless of sample number
#HT 01 and 05 were day 1, Ht 02 and 06 were day 15, hashtag 03 and 07 were day 90, and hashtag 04 and 08 were day 180
mergedPosObjs <- mergedPosObjs %>% filter(!(orig.ident == "COV_09_13_2023_P1" & MULTI_ID %in% c("HTO-0255", "HTO-0256"))) #remove test leukopaks that were sequenced with this donor
mergedPosObjs <- mergedPosObjs%>%mutate(Timepoint = case_when(MULTI_ID %in% c("HTO-0251","HTO-0255") ~ "Day 0",
                                                       MULTI_ID %in% c("HTO-0252", "HTO-0256")  ~ "Day 15",
                                                       MULTI_ID %in% c("HTO-0253", "HTO-0257")  ~ "Day 90",
                                                       MULTI_ID %in% c("HTO-0254", "HTO-0258")  ~ "Day 180",
                                                       TRUE ~ "Error in Hashtag Call"
                                                       
              ))

#add sample names now too
#define which half of the hashtags a donor was included in for a given pool
firstHalf <- c("HTO-0251", "HTO-0252", "HTO-0253", "HTO-0254")
lastHalf <- c("HTO-0255", "HTO-0256", "HTO-0257", "HTO-0258")
 
#add names
mergedPosObjs <- mergedPosObjs %>% mutate(Subject = case_when(orig.ident %in% c("COV_09_13_2023_P1") & MULTI_ID %in% firstHalf ~ "5351564848",
                                                              
                                                              orig.ident %in% c("COV_09_19_2023_P1") & MULTI_ID %in% firstHalf ~ "5457484948",
                                                              orig.ident %in% c("COV_09_19_2023_P1") & MULTI_ID %in% lastHalf ~ "5048574848",
                                                              
                                                              orig.ident %in% c("COV_09_19_2023_P2","COV_09_19_2023_P3") & MULTI_ID %in% firstHalf ~ "4951574848",
                                                              orig.ident %in% c("COV_09_19_2023_P2","COV_09_19_2023_P3") & MULTI_ID %in% lastHalf ~ "5357484948",
                                                              
                                                              
                                                              orig.ident %in% c("COV_09_20_2023_P1") & MULTI_ID %in% firstHalf ~ "4954554848",
                                                              orig.ident %in% c("COV_09_20_2023_P1") & MULTI_ID %in% lastHalf ~ "4957484948",
                                                              
                                                              orig.ident %in% c("COV_09_20_2023_P2","COV_09_20_2023_P3") & MULTI_ID %in% firstHalf ~ "4848544848",
                                                              orig.ident %in% c("COV_09_20_2023_P2","COV_09_20_2023_P3") & MULTI_ID %in% lastHalf ~ "5053564848",
                                                              
                                                              
                                                              orig.ident %in% c("COV_09_21_2023_P1") & MULTI_ID %in% firstHalf ~ "4955534848",
                                                              orig.ident %in% c("COV_09_21_2023_P1") & MULTI_ID %in% lastHalf ~ "4957564848",
                                                              
                                                              orig.ident %in% c("COV_09_21_2023_P2") & MULTI_ID %in% firstHalf ~ "4953494948",
                                                              orig.ident %in% c("COV_09_21_2023_P2") & MULTI_ID %in% lastHalf ~ "5553564848",
                                                              
                                                              
                                                              orig.ident %in% c("COV_09_26_2023_P1") & MULTI_ID %in% firstHalf ~ "5150564848",
                                                              orig.ident %in% c("COV_09_26_2023_P1") & MULTI_ID %in% lastHalf ~ "5557544848",
                                                              
                                                              orig.ident %in% c("COV_09_26_2023_P2") & MULTI_ID %in% firstHalf ~ "5050484948",
                                                              orig.ident %in% c("COV_09_26_2023_P2") & MULTI_ID %in% lastHalf ~ "5048544848",
                                                              
                                                              
                                                              orig.ident %in% c("COV_09_27_2023_P1") & MULTI_ID %in% firstHalf ~ "5750564848",
                                                              orig.ident %in% c("COV_09_27_2023_P1") & MULTI_ID %in% lastHalf ~ "5755544848",
                                                              
                                                              orig.ident %in% c("COV_09_27_2023_P2") & MULTI_ID %in% firstHalf ~ "4950544848",
                                                              orig.ident %in% c("COV_09_27_2023_P2") & MULTI_ID %in% lastHalf ~ "5249544848",
                                                              
                                                              
                                                              orig.ident %in% c("COV_09_28_2023_P1", "COV_09_28_2023_P3") & MULTI_ID %in% firstHalf ~ "5054574848",
                                                              orig.ident %in% c("COV_09_28_2023_P1", "COV_09_28_2023_P3") & MULTI_ID %in% lastHalf ~ "5456544848",
                                                              
                                                              orig.ident %in% c("COV_09_28_2023_P2") & MULTI_ID %in% firstHalf ~ "5356534848",
                                                              orig.ident %in% c("COV_09_28_2023_P2") & MULTI_ID %in% lastHalf ~ "5556484948",
                                                              
                                                              
                                                              orig.ident %in% c("COV_10_05_2023_P1") & MULTI_ID %in% firstHalf ~ "5553554848",
                                                              
                                                              
                                                              orig.ident %in% c("COV_10_06_2023_P1") & MULTI_ID %in% firstHalf ~ "5657574848",
                                                              
                                                              orig.ident %in% c("COV_10_06_2023_P2") & MULTI_ID %in% firstHalf ~ "4856554848",
                                                              orig.ident %in% c("COV_10_06_2023_P2") & MULTI_ID %in% lastHalf ~ "5653544848",
                                                              
                                                              TRUE ~ "Error in labelling donor"
                                                              ))

#Now let's add booster annotation:
o <- c("5048574848",
       "4957484948",
       "4848544848",
       "4957564848",
       "5351564848")

p <- c("5457484948",
       "4951574848",
       "4954554848",
       "4955534848",
       "5557544848")

op <- c("5357484948",
        "5053564848",
        "4953494948",
        "5553564848",
        "5150564848")

oi <- c("5050484948",
        "5750564848",
        "5755544848",
        "5054574848",
        "5456544848")

pi <- c("5048544848",
        "4950544848",
        "5249544848",
        "5356534848",
        "5556484948")

opi <- c("5553554848",
         "5657574848",
         "4856554848",
         "5653544848")

#add the annotation
mergedPosObjs <- mergedPosObjs%>%mutate(Booster = case_when(Subject %in% c(o, oi) ~ "Omicron",
                                                            Subject %in% c(p, pi) ~ "Prototype",
                                                            Subject %in% c(op, opi) ~ "Omicron And Prototype",
                                                            TRUE ~ "Error in annotating booster dose"))
                                                              
mergedPosObjs <- mergedPosObjs%>%mutate(Infection = case_when(Subject %in% c(oi, pi, opi) ~ "Y",
                                                              Subject %in% c(o, p, op) ~ "N",
                                                              TRUE ~ "Error in labelling infection status"
  
))

message("Counts for cells per subject: \n",paste0(capture.output(plyr::count(mergedPosObjs@meta.data$Subject)), collapse = "\n"))
message("Counts for cells per timepoint for all pools: \n",paste0(capture.output(plyr::count(mergedPosObjs@meta.data$Timepoint)), collapse = "\n"))
message("Counts for cells per booster dose: \n",paste0(capture.output(plyr::count(mergedPosObjs@meta.data$Booster)), collapse = "\n"))
message("Counts for cells by infection status: \n",paste0(capture.output(plyr::count(mergedPosObjs@meta.data$Infection)), collapse = "\n"))    

#edit run info file
run_info$neg.drops[run_number] <- dim(mergedNegObjs)[2]
run_info$wt.cells[run_number] <- dim(mergedPosObjs)[2]

saveRDS(mergedPosObjs, file = here::here("04_Analysis","data_objects","01_build_seurat","MergedSeuratObject_p.rds"))
saveRDS(mergedNegObjs, file = here::here("04_Analysis","data_objects","01_build_seurat","MergedSeuratObject_n.rds"))

write.csv(run_info,file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"), row.names = FALSE)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()