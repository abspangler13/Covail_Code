#test out cell selector 

library(Seurat)

A316 <- readRDS("~/Desktop/A316_HA.rds")
#A316 <- readRDS("/Volumes/vrc_vip/Abby/Experiment_316/analysis/code/04_HA/A316_HA.rd")

plot <- DimPlot(object = A316,reduction = "harmony.rna.umap")
# Follow instructions in the terminal to select points
memB <- CellSelector(plot = plot)
PB <- CellSelector(plot = plot)
junk <- CellSelector(plot = plot)
junk2 <- CellSelector(plot = plot)

junk<-c(junk,junk2)

PB <- setdiff(PB,junk)
memB <- setdiff(memB,junk)

cell_selector <- list(memB,PB,junk)
names(cell_selector) <- c("memB","PB","junk")
saveRDS(cell_selector, "cell_selector.rds")
