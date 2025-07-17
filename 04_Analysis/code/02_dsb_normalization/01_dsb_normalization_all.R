library(Seurat)
library(here)
library(tidyverse)
library(sessioninfo)
library(dsb)
library(limma)

# Read in negative and positive objects as determined by hashtag demultiplexting. Load in good object which is after QC.
mergedNegObject <- readRDS(file = here::here("04_Analysis","data_objects","01_build_seurat","MergedSeuratObject_n.rds"))
mergedPosObject <- readRDS(file = here::here("04_Analysis","data_objects","01_build_seurat","MergedSeuratObject_p.rds"))

DefaultAssay(mergedPosObject) <- "Prot"

message("Dimensions of merged positive object before setting upper prot limit: ", dim(mergedPosObject)[1], ", ", dim(mergedPosObject)[2])
# [1] 33538 66204

######QC protein libraries########
prot.mult = (3*mad(mergedPosObject$log_nCount_Prot))
prot.upper = median(mergedPosObject$log_nCount_Prot) + prot.mult
mergedPosObject <- mergedPosObject %>% subset(subset = log_nCount_Prot < prot.upper)

message("Dimensions of merged positive object after filtering out cells with prot counts 3x above MAD: ", dim(mergedPosObject)[1], ", ", dim(mergedPosObject)[2])
# [1]  33538 6506


#umap of before
mergedPosObject <- ScaleData(mergedPosObject, features = rownames(mergedPosObject))
mergedPosObject <- RunPCA(mergedPosObject, assay = "Prot", slot = "data", features = rownames(mergedPosObject), reduction.name = "apca")

mergedPosObject <- RunUMAP(mergedPosObject, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
pdf(file = here::here("04_Analysis","plots","02_dsb_normalization","DimPlot_prot_UMAP_all_before_normalization.pdf"))
DimPlot(mergedPosObject, reduction = "prot.umap", label = TRUE)
DimPlot(mergedPosObject, reduction = "prot.umap", label = TRUE, group.by = "orig.ident")
DimPlot(mergedPosObject, reduction = "prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

############################################Do DSB normalization ########################################
# https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html 

isotypes = rownames(mergedPosObject@assays$Prot@counts)[grep("iso",rownames(mergedPosObject@assays$Prot@counts))]

normalized_dsb_matrix_sm = DSBNormalizeProtein(cell_protein_matrix = as.matrix(mergedPosObject@assays$Prot@counts),
                                               empty_drop_matrix = as.matrix(mergedNegObject@assays$Prot@counts), use.isotype.control = TRUE, isotype.control.name.vec = isotypes)

# #####Run second linear correction using Limma to reduce batch effect######
# normalized_dsb_matrix_sm = removeBatchEffect(x = normalized_dsb_matrix_sm, batch = mergedPosObject$run, batch2 = mergedPosObject$Subject)

#Save assay data into normalized data slot 
mergedPosObject <- SetAssayData(mergedPosObject, slot = "data", new.data = normalized_dsb_matrix_sm, assay = "Prot")

saveRDS(mergedPosObject, file = here::here("04_Analysis","data_objects","02_dsb_normalization","MergedSeuratObject_positive_dsbnormalized.rds"))

# mergedPosObject <-mergedPosObject

#do we need to set variable features here? VariableFeatures(bm) <- rownames(bm[["ADT"]])
mergedPosObject <- ScaleData(mergedPosObject, features = rownames(mergedPosObject))
mergedPosObject <- RunPCA(mergedPosObject, assay = "Prot", slot = "data", features = rownames(mergedPosObject), reduction.name = "apca")

mergedPosObject <- RunUMAP(mergedPosObject, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
pdf(file = here::here("04_Analysis","plots","02_dsb_normalization","DimPlot_prot_UMAP_all_after_normalization.pdf"))
DimPlot(mergedPosObject, reduction = "prot.umap", label = TRUE)
DimPlot(mergedPosObject, reduction = "prot.umap", label = TRUE, group.by = "orig.ident")
DimPlot(mergedPosObject, reduction = "prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
