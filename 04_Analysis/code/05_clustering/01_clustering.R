library(Seurat)
library(sessioninfo)
library(harmony)
library(tidyseurat)
library(tidyverse)
library(ggridges)
library(here)
library(devtools)

#set a seed for reproducibility when running on HPC
set.seed(1)

print("R Version:")
print(R.version)

###################################### Clustering based on Prot and RNA transcripts #################################
sObj <- readRDS(file = here::here("04_Analysis","data_objects","04_probe","BindingDataCorrected_COVAIL_SeuratObj.rds"))

## Cluster on only RBD+ B-cells
#Note from 240228 Rory - Since we only care about prototype and omicron, I've shifted this to only proto+ and/or omi+ cells
rbdPos <- sObj %>% filter(adj.ProtoOmi != "Proto-Omi-")

## Did clustering in integration script, but will re-do it here because we dropped junk cells in 04_HA. Also want to drop HA- cells before clustering
DefaultAssay(rbdPos) <- "RNA"
rbdPos <- NormalizeData(rbdPos, normalization.method = "LogNormalize")
rbdPos <- FindVariableFeatures(rbdPos, selection.method = "vst", nfeatures = 600)
VariableFeatures(rbdPos) <- grep("IG[HKL]V|TRBV", VariableFeatures(rbdPos), invert = TRUE, value = TRUE)
rbdPos <- ScaleData(rbdPos)
rbdPos <- RunPCA(rbdPos, features = VariableFeatures(rbdPos))
rbdPos <- RunHarmony(rbdPos, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")

#make umaps of both PCs and harmony corrected PCs
rbdPos <- RunUMAP(rbdPos, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA") #do I need to change the number of dims based on the elbow plot?
rbdPos <- RunUMAP(rbdPos, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

#we do have demultiplexed subject data, but our sequencing data was done
#by date of sort
pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_rna_UMAP_RBDPositive.pdf"))
VariableFeaturePlot(rbdPos)
ElbowPlot(rbdPos)
DimPlot(rbdPos, reduction = "rna.umap",group.by = "Subject")
DimPlot(rbdPos, reduction = "harmony.rna.umap",group.by = "Subject")
# DimPlot(rbdPos, reduction = "rna.umap",group.by = "run") #we don't have a run date- the run date is basically orig.ident
# DimPlot(rbdPos, reduction = "harmony.rna.umap",group.by = "run")
DimPlot(rbdPos, reduction = "rna.umap", label = TRUE) #consider removing label=TRUE because it gets in the way of the plot
DimPlot(rbdPos, reduction = "rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(rbdPos, reduction = "harmony.rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(rbdPos, reduction = "rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

# Run on Prot individually #
DefaultAssay(rbdPos) <- "Prot"
ProtFeatures <- rownames(rbdPos@assays$Prot@data)
VariableProtFeatures <- grep("Ig[M|G|D|A]", ProtFeatures, invert = TRUE, value = TRUE)
VariableProtFeatures <- ProtFeatures
rbdPos <- ScaleData(rbdPos, features = VariableProtFeatures)
rbdPos <- RunPCA(rbdPos, assay = "Prot", slot = "data", features = rownames(rbdPos), reduction.name = "apca")
rbdPos <- RunHarmony(rbdPos, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

#make umaps of both PCs and harmony corrected PCs - do we need to change the number of dimensions?
rbdPos <- RunUMAP(rbdPos, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
rbdPos <- RunUMAP(rbdPos, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_prot_UMAP.pdf"))
ElbowPlot(rbdPos)
DimPlot(rbdPos, reduction = "prot.umap",group.by = "Subject")
DimPlot(rbdPos, reduction = "harmony.prot.umap",group.by = "Subject")
# DimPlot(rbdPos, reduction = "prot.umap",group.by = "run")
# DimPlot(rbdPos, reduction = "harmony.prot.umap",group.by = "run")
DimPlot(rbdPos, reduction = "prot.umap", group.by = "orig.ident")
DimPlot(rbdPos, reduction = "harmony.prot.umap", group.by = "orig.ident")
DimPlot(rbdPos, reduction = "harmony.prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

# Do multi-modal clustering
rbdPos <- FindMultiModalNeighbors(rbdPos, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")
rbdPos <- FindClusters(rbdPos, graph.name = "harmony.snn", algorithm = 3, resolution = 0.4, verbose = FALSE)
rbdPos <- RunUMAP(rbdPos, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP.pdf"))
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE) + theme(aspect.ratio = 1)
# DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "broad.celltype", label = TRUE)
# DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "run", split.by = "run", ncol = 2)
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "orig.ident")
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "Subject")
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering.pdf"))
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Subject", ncol = 4)
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "orig.ident", ncol = 2)
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Booster", ncol = 4)
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Infection", ncol = 4)
dev.off()

table(rbdPos$harmony.snn_res.0.4)
# 0    1   10   11   12   13    2    3    4    5    6    7    8    9 
# 8940 5495  146  133   99   35 3111 2973 2920 1374 1314  875  749  160 

saveRDS(rbdPos, file = here::here("04_Analysis","data_objects","05_clustering","COVAIL_ClusteredSeuratObject.rds"))

# look at diff exp protein markers all clusters with res 0.4 #
Prot <- data.frame(t(as.matrix(rbdPos@assays$Prot@data)))
Prot$CELL <- rownames(Prot)
head(Prot)
meta <- rbdPos@meta.data[, c("CELL", "harmony.snn_res.0.4", "Population", "Timepoint", "Subject", "c_call")]
Prot.meta <- left_join(Prot, meta, by = "CELL")

############# DSb analysis PDF #############
dsb.analysis <- Prot.meta[,c(1:59)] %>% filter(!(harmony.snn_res.0.4 %in% c("13")))
Proteins <- colnames(dsb.analysis)[1:53]

dsb.plot.list <- list()

pdf(file = here::here("04_Analysis","plots","05_clustering","Protein_ridge_plots_before.pdf"))
for ( i in 1:length(Proteins)) {
  
  dsb.plot.list[[i]] <- dsb.analysis %>%
    ggplot(aes_string(x = Proteins[i], y = "harmony.snn_res.0.4")) + 
    geom_density_ridges(alpha = 0.6) + 
    theme_classic()
  
  plot(dsb.plot.list[[i]])
  
}
dev.off()

#below we remove cluster 11 from our analysis, but it's just mildly suspicious- let's do DE
markers <- FindAllMarkers(rbdPos, assay = "RNA") %>%
  filter(abs(avg_log2FC) > 0.5, !grepl("IGH", rownames(.)))

write.csv(markers, here::here("04_Analysis", "data_objects", "05_clustering", "PreFiltering_DEAnalysis.csv"))

#make a probe plot on top of the dimplot
pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering_probe_label.pdf"))
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "adj.ProtoOmi")
DimPlot(rbdPos, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
dev.off()

#notes on some of the ugly clusters (compared to Locus, clusters 10 and 11 have swapped names- not totally sure why)
#10 has high CD 14 (likely macrophage)
#12 and 7 express CD3
#7 is low for CD19
#11 seems to have a split BAFFR+ and BAFFR- population? and also CD38 high- keep for now
#4 is a bit high for CD14 and lower for CD19 and CD20
#setting a cutoff for CD3 and CD4 that might remove some of the cluster 10 ickiness
print("Removing the following clusters: 4, 7, 10, 12")
rbdPos.G <- rbdPos %>% filter(!(harmony.snn_res.0.4 %in% c("4" ,"7", "10", "12"))) #running locally
#removed.tab <- rbdPos %>% filter(harmony.snn_res.0.4 %in% c("7", "9", "10", "8", )) %>%

## Re-do clustering with just memory cells
DefaultAssay(rbdPos.G) <- "RNA"
rbdPos.G <- NormalizeData(rbdPos.G, normalization.method = "LogNormalize")
rbdPos.G <- FindVariableFeatures(rbdPos.G, selection.method = "vst", nfeatures = 600)
VariableFeatures(rbdPos.G) <- grep("IG[HKL]V|TRBV", VariableFeatures(rbdPos.G), invert = TRUE, value = TRUE)
rbdPos.G <- ScaleData(rbdPos.G)
rbdPos.G <- RunPCA(rbdPos.G, features = VariableFeatures(rbdPos.G))
rbdPos.G <- RunHarmony(rbdPos.G, group.by.vars = "orig.ident", reduction = "pca", reduction.save = "harmony")
rbdPos.G <- RunUMAP(rbdPos.G, reduction = "pca", reduction.name = "rna.umap", dims = 1:30, assay = "RNA")
rbdPos.G <- RunUMAP(rbdPos.G, reduction = "harmony", reduction.name = "harmony.rna.umap", dims = 1:30, assay = "RNA")

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_rna_UMAP_G.pdf"))
VariableFeaturePlot(rbdPos.G)
ElbowPlot(rbdPos.G)
DimPlot(rbdPos.G, reduction = "rna.umap",group.by = "Subject")
DimPlot(rbdPos.G, reduction = "harmony.rna.umap",group.by = "Subject")
# DimPlot(rbdPos.G, reduction = "rna.umap",group.by = "run")
# DimPlot(rbdPos.G, reduction = "harmony.rna.umap",group.by = "run")
DimPlot(rbdPos.G, reduction = "rna.umap", label = TRUE)
DimPlot(rbdPos.G, reduction = "rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(rbdPos.G, reduction = "harmony.rna.umap", label = TRUE, group.by = "orig.ident")
DimPlot(rbdPos.G, reduction = "rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
DimPlot(rbdPos.G, reduction = "harmony.rna.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

DefaultAssay(rbdPos.G) <- "Prot"
ProtFeatures <- rownames(rbdPos.G@assays$Prot@data)
#VariableProtFeatures <- ProtFeatures[- c(47, 48, 53, 56, 58)]
VariableProtFeatures <- grep("Ig[M|G|D|A]", ProtFeatures, invert = TRUE, value = TRUE)
rbdPos.G <- ScaleData(rbdPos.G, features = VariableProtFeatures)
rbdPos.G <- RunPCA(rbdPos.G, assay = "Prot", slot = "data", features = rownames(rbdPos.G), reduction.name = "apca")
rbdPos.G <- RunHarmony(rbdPos.G, group.by.vars = "orig.ident", reduction = "apca", reduction.save = "harmony.prot")

#make umaps of both PCs and harmony corrected PCs
rbdPos.G <- RunUMAP(rbdPos.G, reduction = "apca", dims = 1:18, assay = "Prot", reduction.name = "prot.umap", reduction.key = "protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)
rbdPos.G <- RunUMAP(rbdPos.G, reduction = "harmony.prot", dims = 1:18, assay = "Prot", reduction.name = "harmony.prot.umap", reduction.key = "har_protUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_prot_UMAP_memB_G.pdf"))
ElbowPlot(rbdPos.G)
DimPlot(rbdPos.G, reduction = "prot.umap",group.by = "Subject")
DimPlot(rbdPos.G, reduction = "harmony.prot.umap",group.by = "Subject")
# DimPlot(rbdPos.G, reduction = "prot.umap",group.by = "run")
# DimPlot(rbdPos.G, reduction = "harmony.prot.umap",group.by = "run")
DimPlot(rbdPos.G, reduction = "prot.umap", group.by = "orig.ident")
DimPlot(rbdPos.G, reduction = "harmony.prot.umap", group.by = "orig.ident")
DimPlot(rbdPos.G, reduction = "harmony.prot.umap", ncol=4, group.by = "Subject", split.by = "Subject")
dev.off()

# Combine Prot and RNA clusters #
rbdPos.G <- FindMultiModalNeighbors(rbdPos.G, reduction.list = list("harmony", "harmony.prot"), dims.list = list(1:15, 1:18), modality.weight.name = "harmony.weight", snn.graph.name = "harmony.snn", knn.graph.name = "harmony.knn", weighted.nn.name = "harmony.weighted.wnn")
rbdPos.G <- FindClusters(rbdPos.G, graph.name = "harmony.snn", algorithm = 3, resolution = 0.4, verbose = FALSE)

rbdPos.G <- RunUMAP(rbdPos.G, nn.name = "harmony.weighted.wnn", reduction.name = "harmony.wnn.umap", reduction.key = "harmony.wnnUMAP_", n.neighbors = 40, min.dist = 0.3, local.connectivity = 3, spread = 3)

table(rbdPos.G$harmony.snn_res.0.4)

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_G.pdf"))
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "orig.ident")
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "Subject")
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "Subject", split.by = "Subject", ncol = 4)
dev.off()

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering_G.pdf"))
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE) + theme(aspect.ratio = 1)
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", label = TRUE)
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "Subject", ncol = 4)
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4", split.by = "orig.ident", ncol = 2)
dev.off()

#Dimplot specifically by timepoint to see how populations change over time
# rbdPos.G$Timepoint <- factor(rbdPos.G$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))
# rbdPos.G$seurat_clusters <- factor(rbdPos.G$seurat_clusters, levels = c(unique(rbdPos.G$seurat_clusters)))

#Seurat is being weird about the legend, so I'll edit it myself by readding cell embeddings
metadata <- rbdPos.G@meta.data
metadata <- cbind(metadata, rbdPos.G@reductions$harmony.wnn.umap@cell.embeddings)
metadata$Timepoint <- factor(metadata$Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180"))

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering_G_OverTime.pdf"), width = 15, height=4)
ggplot(metadata[metadata$Infection == "N",], aes(x=harmonywnnUMAP_1, y=harmonywnnUMAP_2))+
  geom_point(aes(color= seurat_clusters, fill=seurat_clusters), size=0.3)+
  ggtitle("UMAP Clustering Over Time")+
  facet_grid(cols = vars(Timepoint))+
  theme_classic()+
  guides(color = guide_legend(override.aes = list(size=3)))
dev.off()

#####
saveRDS(rbdPos.G, file = here::here("04_Analysis","data_objects","05_clustering","COVAIL_ClusteredSeuObj_CleanedAndReclustered.rds"))

#create an after version of the protein plots
# look at diff exp protein markers all clusters with res 0.4 #
Prot <- data.frame(t(as.matrix(rbdPos.G@assays$Prot@data)))
Prot$CELL <- rownames(Prot)
head(Prot)
meta <- rbdPos.G@meta.data[, c("CELL", "harmony.snn_res.0.4", "Population", "Timepoint", "Subject", "c_call")]
Prot.meta <- left_join(Prot, meta, by = "CELL")

pdf(file = here::here("04_Analysis","plots","05_clustering","DimPlot_UMAP_clustering_probe_label_G.pdf"))
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "adj.ProtoOmi")
DimPlot(rbdPos.G, reduction = "harmony.wnn.umap", group.by = "harmony.snn_res.0.4")
dev.off()

############# DSb analysis PDF #############
dsb.analysis <- Prot.meta[,c(1:59)] %>% filter(!(harmony.snn_res.0.4 %in% c("13")))
Proteins <- colnames(dsb.analysis)[1:53]

dsb.plot.list <- list()

pdf(file = here::here("04_Analysis","plots","05_clustering","Protein_ridge_plots_aftercleaning.pdf"))
for ( i in 1:length(Proteins)) {

  dsb.plot.list[[i]] <- dsb.analysis %>%
    ggplot(aes_string(x = Proteins[i], y = "harmony.snn_res.0.4")) +
    geom_density_ridges(alpha = 0.6) +
    theme_classic()

  plot(dsb.plot.list[[i]])

}
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()