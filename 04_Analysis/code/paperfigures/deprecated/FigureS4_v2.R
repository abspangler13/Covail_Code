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

#scale the data
seuObjNaive <- ScaleData(seuObjNaive, features <- rownames(seuObjNaive))

#do DE analysis at a single cell level
markers <- FindAllMarkers(seuObjNaive, assay = "RNA", group.by = "ClusterLabel") %>%
  filter(abs(avg_log2FC) > 1, p_val_adj < 0.05, !grepl("IGH", rownames(.)))

write_xlsx(markers, here::here("04_Analysis", "data_objects", "paperfigures", "Figure S4", "DEAnalysis_SignificantMarkers.xlsx"))

#read in csv if i want to work locally
# markers <- read_xlsx(here::here("04_Analysis", "data_objects", "paperfigures", "Figure S4", "DEAnalysis_SignificantMarkers.xlsx"))

#whittle down list of genes
markersList <-  markers %>% group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 20) %>% slice_max(order_by = avg_log2FC, n = 20)

#####pseudobulk the data to make it plottable
#AggregateExpression in Seurat is acting weird. The resulting values don't always match the patterns we expect to see from the data using VlnPlot()
#I think cluster size is biasing results so that genes in large clusters that aren't proportionally very upregulated seem much more so
#let's manually pseudobulk by calculating mean expression per cluster
normalizedRNA <- as.data.frame(t(seuObjNaive@assays$RNA@scale.data)) %>% tibble::rownames_to_column() %>% rename("CELL" = "rowname") %>%
                  mutate(ClusterLabel = seuObjNaive$ClusterLabel[match(CELL, seuObjNaive$CELL)]) %>% select(!CELL) %>%
                  group_by(ClusterLabel) %>% summarise_all(mean)
rownames(normalizedRNA) <- normalizedRNA$ClusterLabel
pseudoObj <- t(normalizedRNA)[-1,]


#normalize the data myself
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

class(pseudoObj) <- "numeric"

#normalize data
pseudoObj <- apply(pseudoObj, 1, cal_z_score)

####let's choose the most significant genes from the list to include in the data- maybe 10?
#plotMarkers <- markers %>% filter(avg_log2FC > 0) %>% group_by(cluster) %>% arrange(p_val_adj, desc(avg_log2FC)) %>% slice(1:20)
pseudoObj <- pseudoObj[,colnames(pseudoObj) %in% markersList$gene]

#pseudoOrder <- pseudoObj[c("Atypical","Acute Activated", "Intermediate", "Resting IgG", "Resting IgA", "Plasmablast-like", "Naive"),]
pseudoOrder <- pseudoObj
pseudoOrder <- pseudoOrder[c("Plasmablast-like", "Naive", "Resting IgA", "Resting IgG", "Atypical", "Intermediate", "Acute Activated"),]

####create two annotations- one for most upregulated gene and one for second
gene <- markers %>% filter(gene %in% markersList$gene) %>%
  group_by(gene) %>% arrange(desc(avg_log2FC), p_val_adj) %>%
  mutate(ann = cluster[1], ann2 = cluster[2]) %>% slice(1)%>%
  select(gene, ann, ann2) %>% tibble::column_to_rownames("gene") %>%
  mutate(ann = ifelse(is.na(ann), "None", ann),
         ann2 = ifelse(is.na(ann2), "None", ann2))

####pheatmaps plotting
class(pseudoOrder) <- "numeric"
svg(filename = here::here("04_Analysis", "plots", "paperfigures", "Figure S3", "DEGHeatmap_method2.svg"), height = 13, width = 4.5)
pheatmap(pseudoOrder, cluster_rows = FALSE,
         annotation_col = gene,
         annotation_colors = shortColors,
         border_color = "gray45",
         fontsize =6.5
)
dev.off()