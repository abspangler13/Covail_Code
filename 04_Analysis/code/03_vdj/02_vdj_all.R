library(tidyverse)
library(Seurat)
library(tidyseurat)
library(sessioninfo)
library(here)

####script to add vdj info to seurat object 
run_info <- read.csv(file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"))
runNum <- 1

#load seurat object 
dsbObj <- readRDS(file = here::here("04_Analysis","data_objects","02_dsb_normalization","MergedSeuratObject_positive_dsbnormalized_BCellAnnotated.rds"))

#load final product of vdj
vdj <- read.csv(file = here::here("04_Analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))

# remove duplicated columns #
vdj <- vdj %>% distinct(CELL,.keep_all=TRUE)
rownames(vdj) <- vdj$CELL

# see number of cells with no Pr or no HT #
apply(dsbObj@meta.data, 2, FUN=function(x) length(which(is.na(x))))
apply(vdj, 2, FUN=function(x) length(which(is.na(x))))

#remove unwanted columns from vdj meta data
vdj <- vdj[,-which(colnames(vdj)%in%c("rev_comp","productive",
                                                        "v_cigar",
                                                        "d_cigar",
                                                        "j_cigar",
                                                        "np1_length",
                                                        "np2_length",
                                                        "v_sequence_start",
                                                        "v_sequence_end",
                                                        "v_germline_start",
                                                        "v_germline_end",
                                                        "d_sequence_start",
                                                        "d_sequence_end",
                                                        "d_germline_start",
                                                        "d_germline_end",
                                                        "j_sequence_start",
                                                        "j_sequence_end",
                                                        "j_germline_start",
                                                        "j_germline_end",
                                                        "chain",
                                                        "v_gene",
                                                        "d_gene",
                                                        "j_gene",
                                                        "c_gene",
                                                        "fwr1",
                                                        "fwr1_nt",
                                                        "cdr1_nt",
                                                        "fwr2",
                                                        "fwr2_nt",
                                                        "cdr2_nt",
                                                        "fwr3",
                                                        "fwr3_nt",
                                                        "cdr3",
                                                        "cdr3_nt",
                                                        "fwr4",
                                                        "fwr4_nt",
                                                        "exact_subclonotype_id",
                                                        "length.l",
                                                        "chain.l",
                                                        "v_gene.l",
                                                        "d_gene.l",
                                                        "j_gene.l",
                                                        "c_gene.l",
                                                        "fwr1.l",
                                                        "fwr1_nt.l",
                                                        "cdr1.l",
                                                        "cdr1_nt.l",
                                                        "fwr2.l",
                                                        "fwr2_nt.l",
                                                        "cdr2.l",
                                                        "cdr2_nt.l",
                                                        "fwr3.l",
                                                        "fwr3_nt.l",
                                                        "cdr3.l",
                                                        "cdr3_nt.l",
                                                        "fwr4.l",
                                                        "fwr4_nt.l",
                                                        "exact_subclonotype_id.l",
                                                        "LC_v_cigar",
                                                        "LC_j_cigar",
                                                        "LC_np1_length",
                                                        "LC_v_sequence_start",
                                                        "LC_v_sequence_end",
                                                        "LC_v_germline_end",
                                                        "LC_j_sequence_start",
                                                        "LC_j_sequence_end",
                                                        "LC_j_germline_start",
                                                        "LC_j_germline_end"))]



run_info$dsbObj[runNum] <- dim(dsbObj)[2]

run_info$all.vdj[runNum] <- dim(vdj)[1]


#add vdj meta data into seurat object 
dsbObj@meta.data$barcode <- rownames(dsbObj@meta.data)
colnames(dsbObj@meta.data)[grep("barcode",colnames(dsbObj@meta.data))] <- "CELL"
run_info$wt.vdj[runNum] <- length(intersect(dsbObj@meta.data$CELL,vdj$CELL))
dsbObj <-left_join(dsbObj,vdj, by="CELL")

apply(dsbObj@meta.data, 2, FUN=function(x) length(which(is.na(x))))

#save whole seurat object
saveRDS(dsbObj, file=here::here("04_Analysis","data_objects","03_vdj","VDJPlusDSBNSeuratObj_final_all.rds"))

#remove cells that don't have vdj and save it as a separate object.
dsbObj.vdj <- dsbObj %>% filter(!is.na(sequence_id))

#concatenate clone_id and subject because some clone_ids are repeated for different subjects
dsbObj.vdj <- dsbObj.vdj %>% mutate(clone_subject_id = paste(Subject, clone_id,sep = '_'))
print(head(dsbObj.vdj$clone_subject_id))

saveRDS(dsbObj.vdj, file=here::here("04_Analysis","data_objects","03_vdj","VDJPlusDSBNSeuratObj_NonVDJRemoved_CloneIDSubjPasted_final_all.rds"))
#dsbObj.vdj <- readRDS(file=here::here("analysis","data_objects","03_vdj","dsbObj_final_vdj_all.rds"))

write.csv(run_info,file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"), row.names=FALSE)

# Make csv with HA and VDJ data for Sarah
dat <- dsbObj.vdj %>% join_features(features=rownames(dsbObj.vdj@assays$Probes)) %>% 
    pivot_wider(names_from=.feature,values_from=.abundance_Probes)
dat <- as.data.frame(dat)

write.csv(dat,file = here::here("04_Analysis","data_objects","03_vdj","VDJ_Probe_metadata.csv"), row.names=FALSE)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


