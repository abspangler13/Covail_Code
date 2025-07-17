#This workflow is based mostly off of Seurat's DESeq2 pseudobulk vignette: https://satijalab.org/seurat/articles/de_vignette

#load in the dependencies
library(Seurat)
library(tidyr)
library(presto)
library(dplyr)
library(tidyseurat)
library(stringr)
library(pheatmap)
library(writexl)
library(readxl)

#load in the data
seuObjNaive <- readRDS(file = here::here("04_Analysis", "data_objects", "06_additional_demultiplexing", "covObj_clustered_demultiplexed.rds")) %>%
  mutate(OfficialBooster = case_when(Booster == "Omicron" ~ "Omicron BA.1 mRNA",
                                     Booster == "Omicron And Prototype" ~ "Prototype + Omicron BA.1 mRNA",
                                     Booster == "Prototype" ~ "Prototype mRNA"),
         ClusterLabel = factor(ClusterLabel, levels = c("Atypical", "Acute Activated", "Intermediate", "Resting IgG", "Resting IgA", "Plasmablast-like", "Naive")),
         Booster = str_replace(Booster, "Omicron", "BA.1"),
         OfficialBooster = factor(OfficialBooster, levels = c("Prototype mRNA", "Prototype + Omicron BA.1 mRNA", "Omicron BA.1 mRNA")),
         Booster = factor(Booster, levels = c("Prototype", "BA.1 And Prototype", "BA.1")))

shortColors <- list(ann = c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                            #"Acute Activated" = "#F46D43",
                            "Acute Activated" = "#f08665",
                            "Intermediate" = "#E6F598",
                            "Resting IgG" = "limegreen",
                            "Resting IgA" = "#3288BD",
                            "Plasmablast-like" = "#6f2da8",
                            "Naive" = "white",
                            "None" = "gray"),
                    ann2 = c("Atypical" = "#D53E4F", #based on RColorBrewer Spectral Palette
                             #"Acute Activated" = "#F46D43",
                             "Acute Activated" = "#f08665",
                             "Intermediate" = "#E6F598",
                             "Resting IgG" = "limegreen",
                             "Resting IgA" = "#3288BD",
                             "Plasmablast-like" = "#6f2da8",
                             "Naive" = "white",
                             "None" = 'gray'))

#scale and normalize data
DefaultAssay(seuObjNaive) <- "RNA"
all.genes <- rownames(seuObjNaive)
seuObjNaive <- ScaleData(seuObjNaive, features = all.genes)

#do DE analysis at a single cell level
markers <- FindAllMarkers(seuObjNaive, assay = "RNA", group.by = "ClusterLabel", min.pct = 0.05) %>%
  filter(abs(avg_log2FC) > 1, p_val_adj < 0.05, !grepl("IGH", rownames(.)))

write_xlsx(markers, here::here("04_Analysis", "data_objects", "paperfigures", "Figure S3", "DEAnalysis_SignificantMarkers.xlsx"))

#read in csv if i want to work locally
# markers <- read_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S3", "DEAnalysis_SignificantMarkers.xlsx"))

#####pseudobulk the data to make it plottable
#pseudoObj <- as.matrix(AggregateExpression(seuObjNaive, assays = "RNA", group.by = c("ClusterLabel"))[["RNA"]])
pseudoObj <- as.matrix(AverageExpression(seuObjNaive, assays = "RNA", group.by = c("ClusterLabel"))[["RNA"]])

####let's choose the most significant genes from the list to include in the data- maybe 10?
plotMarkers <- markers %>% filter(avg_log2FC > 1) %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice(1:20)
pseudoObj <- pseudoObj[rownames(pseudoObj) %in% plotMarkers$gene,]

#function borrowed from https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

pseudoObj <- apply(pseudoObj, 1, cal_z_score)
#pseudoOrder <- pseudoObj[c("Atypical","Acute Activated", "Intermediate", "Resting IgG", "Resting IgA", "Plasmablast-like", "Naive"),]
pseudoOrder <- t(pseudoObj)
pseudoOrder <- pseudoOrder[,c("Atypical", "Acute Activated", "Intermediate","Resting IgG", "Resting IgA", "Naive", "Plasmablast-like")]

####create two annotations- one for most upregulated gene and one for second
gene <- markers %>% filter(gene %in% plotMarkers$gene) %>% filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  group_by(gene) %>% arrange(desc(avg_log2FC)) %>%
  mutate(ann = cluster[1], ann2 = cluster[2]) %>% slice(1)%>%
  select(gene, ann, ann2) %>% tibble::column_to_rownames("gene") %>%
  mutate(ann = ifelse(is.na(ann), "None", ann),
         ann2 = ifelse(is.na(ann2), "None", ann2),
         ann2 = factor(ann, levels = c("Atypical", "Acute Activated", "Intermediate", "Resting IgG", "Resting IgA", "Plasmablast-like", "Naive", "None"),))

geneList <- rownames( gene[with(gene, order(ann, ann2)), ]   )

pseudoOrder2 <- pseudoOrder[order(match(rownames(pseudoOrder), geneList)),drop=FALSE,]

####pheatmaps plotting
svg(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S3", "DEGHeatmap.svg"), height = 13, width = 4.5)
pheatmap(pseudoOrder2, cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_row = gene,
         annotation_colors = shortColors,
         border_color = "gray45",
         fontsize =6.5
)
dev.off()
