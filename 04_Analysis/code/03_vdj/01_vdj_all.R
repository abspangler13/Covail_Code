#I cant get anything with shazam or alakazam to work on Skyline- I'm having version control issues that
#I won't be able to resolve

library(tidyverse)
library(Seurat)
library(readr)
library(plyr)
# library(alakazam)
# library(shazam)
library(devtools)
library(sessioninfo)
library(vroom)
library(here)

# helpful function combined <- AddMetaData(combined, metadata = all_cln_collapsed)
runNum <- 1
run_info <- read.csv(file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"))

####### Add VDJ information ###############
#Read in all annotated vdj from 10X and combine them
vdjData <- vroom(paste0(here::here("03_Immcantation"), "/",list.files(path = here::here("03_Immcantation"),pattern = "*annotations.csv")))[,-1]
#[1] 110488     31

#using VDJpair function to select only barcodes with 1 HC paired with 1 LC#
files <- vdjData
files <- list(files)
temp.file <- list()
paired.files <- list()
i = 1

  temp.file <- files[[i]] %>% dplyr::select(-is_cell, -contig_id, -high_confidence, -full_length, -productive, -raw_clonotype_id, -raw_consensus_id)
  
  temp.file.h <- temp.file %>% dplyr::filter(chain == "IGH")
  #temp.file.k <- temp.file %>% dplyr::filter(chain == "IGK")
  temp.file.l <- temp.file %>% dplyr::filter(chain %in% c("IGL","IGK"))

  #remove duplicate barcodes
  temp.file.h <- temp.file.h[!duplicated(temp.file.h$barcode),]
  #temp.file.k <- temp.file.k[!duplicated(temp.file.k$barcode),]
  temp.file.l <- temp.file.l[!duplicated(temp.file.l$barcode),]

  #temp.HK <- dplyr::left_join(temp.file.h, temp.file.k, by = "barcode")
  temp.HL <- dplyr::left_join(temp.file.h, temp.file.l, by = "barcode")
  
  my.colnames <- c("barcode","length", "chain","v_gene","d_gene",
  "j_gene","c_gene","fwr1","fwr1_nt","cdr1","cdr1_nt","fwr2",
  "fwr2_nt","cdr2","cdr2_nt","fwr3","fwr3_nt","cdr3",
   "cdr3_nt","fwr4","fwr4_nt","reads" ,"umis", "exact_subclonotype_id")
  #colnames(temp.HK) <- c(my.colnames, paste(my.colnames[-1], ".l", sep = ""))
  colnames(temp.HL) <- c(my.colnames, paste(my.colnames[-1], ".l", sep = ""))
  
  #temp.merge <- rbind(temp.HK, temp.HL)
  temp.HL <- temp.HL %>% distinct(barcode, .keep_all = TRUE)

  temp.merge <- temp.HL
  #temp.merge <- temp.merge[!(temp.merge$barcode %in% temp.merge[duplicated(temp.merge$barcode),]$barcode),]
  idents = NULL
  if(length(idents) == length(files)) {
    temp.merge$orig <- paste(idents[i])
  } else {
    temp.merge$orig <- paste("file", i, sep = "_")
  }
  
  paired.files[[i]] <- temp.merge

rtrn <- (do.call(rbind, paired.files))
All.M.vdj <- rtrn

#run_info[nrow(run_info)+1,] <- NA
run_info$vdj.ann.pair[runNum] <- dim(All.M.vdj)[1]

All.M.vdj$CELL <- All.M.vdj$barcode

write.csv(All.M.vdj, file = here::here("04_Analysis","data_objects","03_vdj", "All.M.vdj_PairedSeqs.csv")) 

# Read in clones determined by Immcantation. Used heavy and light chain to determine clones. 
clones <- vroom(paste0(here::here("03_Immcantation"), "/",list.files(path = here::here("03_Immcantation"),pattern = "*heavy_light_clone-pass_germ-pass.tsv")))

colnames(clones)[grep("cell_id",colnames(clones))] <- "CELL"
run_info$VDJCellsByCloneTSV[runNum] <- dim(clones)[1]

#Merging 10X annotation sheet with clonal assignment# left join so only keep cells that have a clone id from immcantation. get rid of remaining contigs from 10x data
All.vdj.M.clones <- left_join(clones, All.M.vdj, by = "CELL") #each naive b-cell_id has unique vdj. daughter cells are clones

# Add in LC info from Immcantation #
#have to put in the output of light chain parse Select function
All.LC <- vroom(paste0(here::here("03_Immcantation"),"/",list.files(path = here::here("03_Immcantation"),pattern = "*light_parse-select.tsv")))

# #addition from Rory on 240216: add in mutations from light chain
# All.LC <- observedMutations(All.LC, sequenceColumn = "sequence_alignment", 
#                             germlineColumn = "germline_alignment", #I'll need to think about how we can run CreateGermlines on LC data- can we do it on the light_parse files?
#                             regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = FALSE)

#All L_ to all column names to distinquish from HC#
colnames(All.LC) <- paste0("LC_", colnames(All.LC))

#Add barcode column
colnames(All.LC)[grep("LC_cell_id",colnames(All.LC))] <- "CELL"

# Remove LCs where more than one per barcode (cell) #
All.LC.singlet <- All.LC[!duplicated(All.LC$CELL),]
run_info$LC.singlet[runNum] <- dim(All.LC.singlet)[1]

#Join LC data with all other data #
All.vdj.M <- left_join(All.vdj.M.clones, All.LC.singlet, by = "CELL")

#remove unwanted *
All.vdj.M$v_call <- gsub("[*].*$", "", All.vdj.M$v_call)
All.vdj.M$LC_v_call <- gsub("[*].*$", "", All.vdj.M$LC_v_call)

#function from shazam. calculates somatic hyper- mutation of vdj. compares sequence of vdj to germline. somatic mutations in vdj. naive bcell has undergone recombination but not somatic hypermutation
# All.vdj.M <- observedMutations(All.vdj.M, sequenceColumn = "sequence_alignment", 
#                                    germlineColumn = "germline_alignment_d_mask", 
#                                    regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)
# 
# count <- observedMutations(All.vdj.M, sequenceColumn = "sequence_alignment", 
#                                germlineColumn = "germline_alignment_d_mask", 
#                                regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = FALSE)
# 
# All.vdj.M$mu_count <- count$mu_count[match(All.vdj.M$sequence_id, count$sequence_id)]
# 
# rm(count)

All.vdj.M  <- All.vdj.M %>% select(-contains("LC_d"))

#final product of vdj
write.csv(All.vdj.M, file = here::here("04_Analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))

# check for where HC, but no LC #
apply(All.vdj.M, 2, FUN=function(x) length(which(is.na(x))))

run_info$clones[runNum] <- length(unique(All.vdj.M$clone_id))

write.csv(run_info,file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"),row.names = FALSE) 

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
