library(Seurat)
library(sessioninfo)
library(harmony)
library(tidyseurat)
library(tidyverse)
library(ggridges)
library(here)

set.seed(1)

#load data
#We usually use azimuth before doing this, but I'm having issues with getting it installed on Skyline, so I'll need to contact the IT team before doing anything
seuObj <- readRDS(file = here::here("04_Analysis","data_objects","05_clustering","COVAIL_ClusteredSeuObj_CleanedAndReclustered.rds"))

#remove cluster 10- this is junk (weird FCRL4 numbers, low cell population anyways)
#remove cluster 9- it seems plasma cell-like in that it expresses high CD38, has two subpopulations of BAFFR expression
seuObj <- seuObj %>% filter(!seurat_clusters %in% c(9, 10))

#Based on Sarah's suggestions, it might be a good idea to recluster again
#We don't care so much about variability in Igs- need to remove them from Prot clustering
#let's put cluster info into metadata and recluster, removing Igs from prot variables and increasing
#resolution to cluster out AM1 subclusters
#add UMAP embeddings and old cluster labels into data
seuObj$SecondClusteringLabel <- seuObj$seurat_clusters

#seuObj <- AddMetaData(seuObj, as.data.frame(Embeddings(seuObj,reduction = "harmony.wnn.umap")))
embeddings <- as.data.frame(Embeddings(seuObj,reduction = "harmony.wnn.umap"))
seuObj@meta.data$SecondUMAP_1 <- embeddings$harmonywnnUMAP_1[match(rownames(embeddings), rownames(seuObj@meta.data))]
seuObj@meta.data$SecondUMAP_2 <- embeddings$harmonywnnUMAP_2[match(rownames(embeddings), rownames(seuObj@meta.data))]


#looks right!
#ggplot(seuObj@meta.data, aes(x = harmonywnnUMAP_1, y = harmonywnnUMAP_2, color = SecondClusteringLabel))+geom_point()+theme_classic()

#recluster
DefaultAssay(seuObj) <- "RNA"
seuObj <- NormalizeData(seuObj, normalization.method = "LogNormalize")
seuObj <- FindVariableFeatures(seuObj, selection.method = "vst", nfeatures = 600)
VariableFeatures(seuObj) <- grep("IG[HKL]V|TRBV", VariableFeatures(seuObj), invert = TRUE, value = TRUE)
seuObj <- ScaleData(seuObj)
seuObj <- RunPCA(seuObj, features = VariableFeatures(seuObj))
seuObj <- RunHarmony(seuObj, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

#make umaps of both PCs and harmony corrected PCs
seuObj <- RunUMAP(seuObj, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA") #do I need to change the number of dims based on the elbow plot?
seuObj <- RunUMAP(seuObj, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

#dimplot to see if harmony was necessary
pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_rna_UMAP_Positive_final_04res.pdf"))
VariableFeaturePlot(seuObj)
ElbowPlot(seuObj)
DimPlot(seuObj, reduction = "rna.umap",group.by = "Subject")
DimPlot(seuObj, reduction = "harmony.rna.umap",group.by = "Subject")
# DimPlot(seuObj, reduction = "rna.umap",group.by = "run") #we don't have a run date- the run date is basically orig.ident
# DimPlot(seuObj, reduction = "harmony.rna.umap",group.by = "run")
DimPlot(seuObj, reduction = "rna.umap", label = TRUE) #consider removing label=TRUE because it gets in the way of the plot
DimPlot(seuObj, reduction = "rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(seuObj, reduction = "harmony.rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(seuObj, reduction = "rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

# Run on Prot individually #
DefaultAssay(seuObj) <- "Prot"
ProtFeatures <- rownames(seuObj@assays$Prot@data)
#VariableProtFeatures <- grep("Ig[M|G|D|A]", ProtFeatures, invert = TRUE, value = TRUE)
VariableProtFeatures <- ProtFeatures
seuObj <- ScaleData(seuObj, features = VariableProtFeatures)
seuObj <- RunPCA(seuObj, assay = "Prot", slot = "data", features = VariableProtFeatures, reduction.name = "apca")
seuObj <- RunHarmony(seuObj, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

#make umaps of both PCs and harmony corrected PCs - do we need to change the number of dimensions?
seuObj <- RunUMAP(seuObj, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
seuObj <- RunUMAP(seuObj, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_prot_UMAP_final_04res.pdf"))
ElbowPlot(seuObj)
DimPlot(seuObj, reduction = "prot.umap",group.by = "Subject")
DimPlot(seuObj, reduction = "harmony.prot.umap",group.by = "Subject")
# DimPlot(seuObj, reduction = "prot.umap",group.by = "run")
# DimPlot(seuObj, reduction = "harmony.prot.umap",group.by = "run")
DimPlot(seuObj, reduction = "prot.umap", group.by = "orig.ident")
DimPlot(seuObj, reduction = "harmony.prot.umap", group.by = "orig.ident")
DimPlot(seuObj, reduction = "harmony.prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

# Do multi-modal clustering
seuObj <- FindMultiModalNeighbors(seuObj, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")
seuObj <- FindClusters(seuObj, graph.name = "harmony.snn", algorithm = 3, resolution = 0.4, verbose = FALSE) #upped the resolution to see what happens
seuObj <- RunUMAP(seuObj, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

# FeaturePlot(seuObj, reduction = "harmony.wnn.umap", features = c("P-IgA"))

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_UMAP_Multimodal_Final_04res.pdf"))
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE) + theme(aspect.ratio = 1)
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "orig.ident")
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "Subject")
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering_final_04res.pdf"))
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Subject", ncol = 4)
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "orig.ident", ncol = 2)
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Booster", ncol = 4)
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Infection", ncol = 4)
dev.off()

saveRDS(seuObj, file = here::here("04_Analysis","data_objects","05_clustering", "res_04","COVAIL_ReclusteredAzimuth_04res.rds"))

#let's look at some plots to see how things stack up
pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_UMAP_clustering_final_04res.pdf"))
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_UMAP_clustering_final_04res_probedata.pdf"))
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "adj.ProtoOmi", label = TRUE)
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_UMAP_clustering_final_04res_plottedagainstsecondplot.pdf"))
ggplot(seuObj@meta.data, aes(x= SecondUMAP_1, y = SecondUMAP_2))+
  geom_point(aes(color = harmony.snn_res.0.4), size = 0.2)+
  ggtitle("Harmony Resolution 0.4 - All Cells")+
  theme_classic()+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

seuObj@meta.data$Timepoint <- factor(seuObj@meta.data$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","UMAP_clustering_final_overtime_04res.pdf"), height = 4, width = 12)
ggplot(seuObj@meta.data[seuObj@meta.data$Infection == "N",], aes(x= SecondUMAP_1, y = SecondUMAP_2))+
  geom_point(aes(color = harmony.snn_res.0.4), size = 0.2)+
  ggtitle("Harmony Resolution 0.4 - Uninfected Over Time")+
  theme_classic()+
  facet_grid(cols = vars(Timepoint))+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","FeaturePlot_UMAP_clustering_04resolution.pdf"), height = 6, width = 10)
FeaturePlot(seuObj, reduction = "harmony.wnn.umap", features = c("P-IgA", "P-IgM"))
dev.off()

#show ridgeplots
# look at diff exp protein markers all clusters with res 0.4 #
Prot <- data.frame(t(as.matrix(seuObj@assays$Prot@data)))
Prot$CELL <- rownames(Prot)
head(Prot)
meta <- seuObj@meta.data[, c("CELL", "harmony.snn_res.0.4", "Population", "Timepoint", "Subject", "c_call")]
Prot.meta <- left_join(Prot, meta, by = "CELL")

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","DimPlot_UMAP_clustering_probe_label_final_04res.pdf"))
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "adj.ProtoOmi")
DimPlot(seuObj, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
dev.off()

############# DSb analysis PDF #############
dsb.analysis <- Prot.meta[,c(1:59)] %>% filter(!(harmony.snn_res.0.4 %in% c("13")))
Proteins <- colnames(dsb.analysis)[1:53]

dsb.plot.list <- list()

pdf(file = here::here("04_Analysis","plots","05_clustering", "res_04","Protein_ridge_plots_final_04res.pdf"))
for ( i in 1:length(Proteins)) {
  
  dsb.plot.list[[i]] <- dsb.analysis %>%
    ggplot(aes_string(x = Proteins[i], y = "harmony.snn_res.0.4")) +
    geom_density_ridges(alpha = 0.6) +
    theme_classic()
  
  plot(dsb.plot.list[[i]])
  
}
dev.off()

#de analysis
#below we remove cluster 11 from our analysis, but it's just mildly suspicious- let's do DE
markers <- FindAllMarkers(seuObj, assay = "RNA") %>%
  filter(abs(avg_log2FC) > 0.5, !grepl("IGH", rownames(.)))

write.csv(markers, here::here("04_Analysis", "data_objects", "05_clustering", "res_04", "FinalDataset_DEAnalysis_04res.csv"))

