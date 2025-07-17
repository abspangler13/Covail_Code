#try azimuth https://satijalab.github.io/azimuth/articles/run_azimuth_tutorial.html

# devtools::install_github("satijalab/seurat-data")
# devtools::install_github("satijalab/azimuth")
# need more than 40G

library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyseurat)
library(tidyverse)
library(sessioninfo)

set.seed(1)

covObj <- readRDS(file = here::here("04_Analysis","data_objects","05_clustering","COVAIL_ClusteredSeuObj_CleanedAndReclustered.rds"))

DefaultAssay(covObj) <- "RNA"
covObj <- RunAzimuth(covObj, reference = "pbmcref", do.adt = TRUE)
saveRDS(covObj,file = here::here("04_Analysis","data_objects","05_clustering","covObj_clustered_azimuth.rds"))

p1 <- DimPlot(covObj, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3, reduction = "harmony.wnn.umap") #not sure if I've been making the clusters wrong, but had to add harmony to wnn.umap
# p2 <- DimPlot(covObj, group.by = "broad.celltype", reduction = "wnn.umap") #we ddidn't sort plasmablasts

pdf(file = here::here("04_Analysis","plots","05_clustering","azimuth_pbmc_predicted_celltype.pdf"))#, width = 36)
p1
dev.off()

print(table(covObj$predicted.celltype.l2))

covObj.sub <- covObj %>% filter(predicted.celltype.l2 %in% c("B intermediate","B memory","B naive"))

print("Successfully subsetted covObj")

p3 <- DimPlot(covObj.sub, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
pdf(file = here::here("04_Analysis","plots","05_clustering","azimuth_pbmc_predicted_celltype_subset.pdf"))#, width = 36)
p3
dev.off()

print("Successfully Dimplotted data subset")

covObj.remove <- covObj %>% filter(!predicted.celltype.l2 %in% c("B intermediate","B memory","B naive"))

print("Successfully subsetted covObj for excluded cells (subsequent dimplot removed)")

#removing the "removed cells" dimplot- the only cells left that aren't intermediate B cells, memory B cells, 
# p5 <- DimPlot(covObj.remove, group.by = "predicted.celltype.l2", label = TRUE, label.size = 3) 
# pdf(file = here::here("04_Analysis","plots","05_clustering","azimuth_pbmc_predicted_celltype_remove.pdf"))#, width = 36)
# p5
# dev.off()
# 
# print("Successfully Dimplotted removed data")

#create table of removed cells 
removed.tab <- covObj.remove %>% join_features(features = rownames(covObj.remove@assays$Probes)) %>% pivot_wider(names_from = .feature, values_from = .abundance_Probes) %>% select(c(1, 10:105))
write.csv(removed.tab,file = here::here("04_Analysis","data_objects","05_clustering","azimuth_removed_cells.csv"))

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()