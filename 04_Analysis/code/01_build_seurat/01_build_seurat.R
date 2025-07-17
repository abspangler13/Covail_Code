library(here)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(sessioninfo)
library(dplyr)
library(tidyr)
library(tidyseurat)
library(scuttle)

set.seed(1)

#I'm having issues with tidyseurat when qsubbing to Locus and have tried reverting it back to 0.7.4 as a potential fix 
message("Tidyseurat version:", packageVersion("tidyseurat"))

#start a working document that contains all cells dropped and what step they were dropped at
removed.tab <- data.frame()

## Define some info for the samples
sample_info <- data.frame(
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
  ),
  cellranger_id = c("VRC-COVAIL-09-13-2023-Pool1",
                    "VRC-COVAIL-09-19-2023-Pool1",
                    "VRC-COVAIL-09-19-2023-Pool2",
                    "VRC-COVAIL-09-19-2023-Pool3",
                    "VRC-COVAIL-09-20-2023-Pool1",
                    "VRC-COVAIL-09-20-2023-Pool2",
                    "VRC-COVAIL-09-20-2023-Pool3",
                    "VRC-COVAIL-09-21-2023-Pool1",
                    "VRC-COVAIL-09-21-2023-Pool2",
                    "VRC-COVAIL-09-26-2023-Pool1",
                    "VRC-COVAIL-09-26-2023-Pool2",
                    "VRC-COVAIL-09-27-2023-Pool1",
                    "VRC-COVAIL-09-27-2023-Pool2",
                    "VRC-COVAIL-09-28-2023-Pool1",
                    "VRC-COVAIL-09-28-2023-Pool2",
                    "VRC-COVAIL-09-28-2023-Pool3",
                    "VRC-COVAIL-10-05-2023-Pool1",
                    "VRC-COVAIL-10-06-2023-Pool1",
                    "VRC-COVAIL-10-06-2023-Pool2"
                    )
)

sample_info$cellranger_raw_dir <- file.path(here::here("02_CellRanger"),
                                        sample_info$cellranger_id,
                                        "outs","raw_feature_bc_matrix")

sample_info$cellranger_filtered_dir <- file.path(here::here("02_CellRanger"),
                                        sample_info$cellranger_id,
                                        "outs","filtered_feature_bc_matrix")

sample_info$protein_csv_file = file.path(here::here("04_Analysis","data_objects","01_build_seurat"),
                             sample_info$sample_id,
                             paste0(sample_info$sample_id,"_sample_specific_oligos.csv"))

sample_info$log.txt = file.path(here::here("04_Analysis","data_objects","01_build_seurat"),
                                sample_info$sample_id,
                                paste0(sample_info$sample_id,"_log.txt"))

for (i in 1:nrow(sample_info)){
  dir.create(file.path(here::here("04_Analysis", "data_objects","01_build_seurat"),
                       sample_info$sample_id[i]))
  dir.create(file.path(here::here("04_Analysis","plots","01_build_seurat"),
                       sample_info$sample_id[i]))
  file.create(file.path(here::here("04_Analysis", "data_objects","01_build_seurat"),
               sample_info$sample_id[i],
               paste0(sample_info$sample_id[i],"_log.txt")))
}

#since we used the same reagents for pretty much every pool (except for the hashtags), we can just copy and paste the same oligo sheet to all sample folders
#basically my feature reference file but with ht nickname and "in.sample" columns
oligoTable <- read.csv(here("04_Analysis","code","01_build_seurat", "my_FeatureReferenceFile_copyforSeuratConstruction.csv"))

#note: Not all samples used every hashtag. Below are hashtags we need to exclude because they don't exist:
#       HT07 and 8 from 9/13
#       HT08 from 9/19, Pool 3
#       HT05 and HT07 from 9/20, Pool 3
#       HT01, HT05, and HT08 from 9/28 pool 3
#       HT05-HT08 for 10/05 Pool 1
#       HT05-HT08 for 10/06 Pool 1

#I previously did this part manually, but that's not very cash money of me, so I'm doing it here instead with a big if-else statement
for(i in 1:nrow(sample_info)){
  placeholder <- oligoTable
  
  if(sample_info$sample_id[i] == "COV_09_13_2023_P1"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT7", "HT8")))
    
  } else if(sample_info$sample_id[i] == "COV_09_19_2023_P3"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT8")))
    
  } else if(sample_info$sample_id[i] == "COV_09_20_2023_P3"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT5", "HT7")))
    
  } else if(sample_info$sample_id[i] == "COV_09_28_2023_P3"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT1", "HT5", "HT8")))
  
  } else if(sample_info$sample_id[i] == "COV_10_05_2023_P1"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT5", "HT6", "HT7", "HT8")))
    
  } else if(sample_info$sample_id[i] == "COV_10_06_2023_P1"){
    placeholder <- placeholder %>% filter(!(ht.key %in% c("HT5", "HT6", "HT7", "HT8")))
    
  } 
  
  write.csv(placeholder, sample_info$protein_csv_file[i], row.names = FALSE)
}
rm(oligoTable)

stopifnot(all(file.exists(sample_info$protein_csv_file)))

QC.plots.b4 <- list()
QC.plots.after <- list()
background.plots <- list()
HTO.vln.plots <- list()
HTO.ridge.plots <- list()
HTO.scatter.plots <- list()
RNA.vln.plots <- list()
RNA.count.vlnplot <- list()

for (i in 1:nrow(sample_info)){
  message("Starting sample ",sample_info$sample_id[i])

   # Read in raw matrix file only so we can get the barcodes that are empty and use them to make a negative seurat object for dsb normalization
  raw.data <- Read10X(data.dir = file.path(sample_info$cellranger_raw_dir[i]))
  message("genes = ", dim(raw.data$`Gene Expression`)[1],", cell barcodes = ",dim(raw.data$`Gene Expression`)[2], " in raw matrix.")
  sample_info$raw.barcodes[i] <- dim(raw.data$`Gene Expression`)[2]

  # Read in cell ranger output. creates a list of two matrices. one is feature barcode 
  # for gene expression and the other is for the CSO custom assay
  filtered.data <- Read10X(data.dir = file.path(sample_info$cellranger_filtered_dir[i]))
  message("genes = ", dim(filtered.data$`Gene Expression`)[1],", cell barcodes = ",dim(filtered.data$`Gene Expression`)[2], " in filtered matrix.")
  sample_info$filtered.barcodes[i] <- dim(filtered.data$`Gene Expression`)[2]
 
  #get cell barcodes for cells and empty drops
  stained_cells = colnames(filtered.data$`Gene Expression`)
  background = setdiff(colnames(raw.data$`Gene Expression`), stained_cells)
  message(length(background), " empty droplets according to CellRanger")
  sample_info$CR.empty.drops[i] <- length(background)
  
  # Create Seurat Object using whole transcriptome gene expression matrix #
  seurat <- CreateSeuratObject(counts = raw.data$`Gene Expression`, project = sample_info$sample_id[i])

  # add indicator for barcodes Cell Ranger called as cells
  seurat@meta.data$CR.drop.class = ifelse(rownames(seurat@meta.data) %in% stained_cells, 'cell', 'background')
  
  # CITE SEQ/Probes and HTOs must be added in separately (as assays) #
  # look st list of oligos to figure out which to include #
  oligoslist <- rownames(as.matrix(raw.data$`Antibody Capture`)) 
  
  #make csv of oligos
  write.csv(oligoslist, file = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"Oligos_list.csv"))

  #load oligo csv that lists which hashtags are present in each sample
  my.csv <- read.csv(file.path(sample_info$protein_csv_file[i]))
  hto.oligos <- my.csv$Oligo[which(my.csv$in.sample ==TRUE & my.csv$assay == "HTO")]
  
  probes.oligos <- my.csv$Oligo[which(my.csv$assay == "Probes")]
  
  prots.oligos <- my.csv$Oligo[which(my.csv$assay == "Prot")]
  
  # Divide up custom assay matrix so we can add hashtags, probes, and proteins as separate assays to the seurat object
  HTO <- as.matrix(raw.data$`Antibody Capture`[hto.oligos,])
  message("HTO assay dim ",dim(HTO)[1]," ",dim(HTO)[2] )
  
  Probes <- as.matrix(raw.data$`Antibody Capture`[probes.oligos,])
  message("Probes assay dim ",dim(Probes)[1]," ",dim(Probes)[2] )
  
  Prot <- as.matrix(raw.data$`Antibody Capture`[prots.oligos,])
  message("Prot assay dim ",dim(Prot)[1]," ",dim(Prot)[2])
  
  #Add matrices as separate assays to seurat object 
  seurat[["HTO"]] <- CreateAssayObject(counts = HTO)
  seurat[["Prot"]] <- CreateAssayObject(counts = Prot)
  seurat[["Probes"]] <- CreateAssayObject(counts = Probes)

  #log normalize for visualization purposes only
  seurat@meta.data$log_nCount_RNA <- log10(Matrix::colSums(seurat@assays$RNA@counts))
  seurat@meta.data$log_nCount_Prot <- log10(Matrix::colSums(seurat@assays$Prot@counts))

  # remove barcodes with no evidence of capture in the experiment. Need to do this so HTO demux runs
  barcodes_keep <- rownames(seurat@meta.data)[which(seurat@meta.data$log_nCount_RNA >0 & seurat@meta.data$log_nCount_Prot > 0)]
  sample_info$zero.capture[i] <- sample_info$raw.barcodes[i] - length(barcodes_keep)

  #save raw object before QC
  seurat_raw<-seurat

  # down select so that QC plots look better 
  # 
  # do another of sliding first then doing drop cells to see if I can recreate the numbers from before
  seurat.old <- seurat %>% arrange(desc(nCount_RNA)) %>% dplyr::slice(1:20000)
  source(file = here::here("04_Analysis","code","QCGenes_FUN.R"))
  seurat.old <- QCGenes(seurat.old)[[1]]
  source(file = here::here("04_Analysis","code","DropCells_FUN.R"))
  seurat_pos.old <- DropCells(seurat.old, min.features = 100, min.hk = 50, max.mito = 15, min.count = 100, max.count = 50000, path.out = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"QC_stats.csv"))
  sample_info$after.old.DropCells[i] <- dim(seurat_pos.old)[2]

  seurat <- subset(seurat, subset = nCount_RNA > 100 & nCount_RNA < 50000)
  message("NumCells after dropping cells with less than 50 and more the 5000 UMIs = ", dim(seurat)[2])
  sample_info$UMI_drop_100_5000[i] <- dim(seurat)[2]

  #Do QC on transcriptome and get rid of bad cells
  #used QCGenes function to calculate percent mito and hk genes. 
  source(file = here::here("04_Analysis","code","QCGenes_FUN.R"))
  seurat<- QCGenes(seurat)[[1]] #adds columns "percent.mito" and "n.exphkgenes" to meta.data

  sample_info$avg.mito[i] <- mean(seurat@meta.data$percent.mito)
  sample_info$avg.hk[i] <- mean(seurat@meta.data$n.exp.hkgenes)
  sample_info$avg.rna[i] <- mean(seurat@meta.data$nCount_RNA)
  sample_info$avg.feature.rna[i] <- mean(seurat@meta.data$nFeature_RNA)
  sample_info$avg.hto[i] <- mean(seurat@meta.data$nCount_HTO)
  sample_info$avg.prot[i] <- mean(seurat@meta.data$nCount_Prot)
  sample_info$avg.probes[i] <- mean(seurat@meta.data$nCount_Probes)
  sample_info$feature.drop[i] <- dim(seurat %>% filter(nFeature_RNA < 100))[2] #having issues with these three lines- 'filter' is not an exported object from namespace:tidyseurat
  sample_info$hk.drop[i] <- dim(seurat %>% filter(n.exp.hkgenes < 50))[2]
  sample_info$mito.drop[i] <- dim(seurat %>% filter(percent.mito > 15))[2]

  #Change this to just use the filter function from tidyseurat. Don't need dropcells function. 
  source(file = here::here("04_Analysis","code","DropCells_FUN.R"))
  seurat_pos <- DropCells(seurat, min.features = 100, min.hk = 50, max.mito = 15, min.count = 100, max.count = 50000, path.out = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],"QC_stats.csv"))
  #test <- seurat %>% filter(nCount_RNA >= 100 & nCount_RNA <= 50000 & percent.mito <= 15 & n.exp.hkgenes >= 50)

  sample_info$after.DropCells[i] <- dim(seurat_pos)[2]
  sample_info$percent.DropCells[i] <- ((dim(seurat)[2]- dim(seurat_pos)[2]) / dim(seurat)[2])*100

  seurat@meta.data$cell_barcode <- colnames(seurat)
  
  #create an object to hold barcodes for cells that were removed
  rnascreenedDroppedCells <- setdiff(paste0(seurat_raw@meta.data$orig.ident, "-", rownames(seurat_raw@meta.data)), paste0(seurat@meta.data$orig.ident, "-", rownames(seurat@meta.data)))
  dropcellsfxnDroppedCells <- setdiff(paste0(seurat@meta.data$orig.ident, "-", rownames(seurat@meta.data)), paste0(seurat_pos@meta.data$orig.ident, "-", rownames(seurat_pos@meta.data)))

  ##use scran too in order to compare
  sce <- seurat@assays$RNA@counts
  is.mito <- grep("MT-",rownames(sce))
  names(is.mito) <- rownames(sce)[is.mito]
  
  df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
  df$cell_barcode <- rownames(df)

 sample_info$med_umi[i] <- summary(df$sum)[3]
  sample_info$mean_umi[i] <- summary(df$sum)[4]
  sample_info$max_umi[i] <- summary(df$sum)[6]
  sample_info$med_feature[i] <- summary(df$detected)[3]
  sample_info$mean_feature[i] <- summary(df$detected)[4]
  sample_info$max_feature[i] <- summary(df$detected)[6]
  sample_info$med_mito[i] <- summary(df$subsets_Mito_percent)[3]
  sample_info$mean_mito[i] <- summary(df$subsets_Mito_percent)[4]
  sample_info$max_mito[i] <- summary(df$subsets_Mito_percent)[6]

  reasons <- perCellQCFilters(df, 
    sub.fields=c("subsets_Mito_percent"))

colSums(as.matrix(reasons))

# > colSums(as.matrix(reasons))
#              low_lib_size            low_n_features high_subsets_Mito_percent 
#                       267                       273                       487 
#                   discard 
#                       693 

  sample_info$low_lib_size[i] <- colSums(as.matrix(reasons))[1]
  sample_info$low_n_features[i] <- colSums(as.matrix(reasons))[2]
  sample_info$high_subsets_Mito_percent[i] <- colSums(as.matrix(reasons))[3]
  sample_info$discard[i] <- colSums(as.matrix(reasons))[4]

  scran.dat <- cbind(df,reasons)
  seurat <- seurat %>% left_join(scran.dat,by = "cell_barcode", copy = TRUE)

  seurat_scran <- seurat %>% filter(discard == "FALSE")
  sample_info$after.ScranDrop[i] <- dim(seurat_scran)[2]
  sample_info$percent.ScranDrop[i] <- ((dim(seurat)[2]- dim(seurat_scran)[2]) / dim(seurat)[2])*100


  #OC on cell barcodes and drop bad cells 
  source(file = here::here("04_Analysis","code","QCplots_FUN.R"))
  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("QCplots_before_DropCells",sample_info$sample_id[i],".pdf")))
  QCplots(seurat)
  dev.off()

  QC.plots.b4[[i]] <- QCplots(seurat)
  
  # make note of which cells are dropped by drop cells and save in raw seurat object
  seurat_raw@meta.data$DropCell = ifelse(rownames(seurat_raw@meta.data) %in% rownames(seurat_pos@meta.data), FALSE, TRUE)
  save(seurat_raw, file=here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0("seurat_raw_",sample_info$sample_id[i],".rds")))

  #look at background drops vs cells
  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("Background_v_cells",sample_info$sample_id[i],".pdf")))
  background.plots[[i]] <- FeatureScatter(seurat_raw,feature1 = "log_nCount_Prot",feature2 = "log_nCount_RNA",group.by = "CR.drop.class") + labs(title = sample_info$sample_id[i])
  print(background.plots[[i]])
  dev.off()

  #QC on background drops. change this to MAD approach from dsb normalization script. 
  seurat_neg <- seurat_raw %>% filter(DropCell==TRUE & CR.drop.class=="background" & log_nCount_Prot>1.5 & log_nCount_Prot <3 & log_nCount_RNA <2.5)
  sample_info$seurat.neg[i] <- dim(seurat_neg)[2]
  
  seurat_neg@meta.data$barcode <- colnames(seurat_neg)
  saveRDS(seurat_neg, file = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_neg.rds")))

  # Normalizing hashtag data by centered log ratio transformation # 
  seurat_pos <- NormalizeData(seurat_pos, assay = "HTO", normalization.method = "CLR")
  
  # Demultiplex HT with MULTIseqDemux For each cell it's designating it with one of the hashtags or as a doublet or negative
  seurat_pos <- MULTIseqDemux(seurat_pos, assay = "HTO", quantile = 0.4)
  
  #not sure if this line was necessary
  Idents(seurat_pos) <- "MULTI_ID"

  # Adding barcodes to meta data as separate column apart from rownames (sommetimes rownames in metadata are converted to numbers) #
  seurat_pos@meta.data$barcode <- colnames(seurat_pos)
  
  #make some plots. Go back and make plos with and with negatives and doublets so we see what dropping them looks like. 
  message("Making some plots")
  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("VlnPlot_nCountHTO_",sample_info$sample_id[i],".pdf")))
  HTO.vln.plots[[i]] <- VlnPlot(seurat_pos, features = "nCount_HTO", group.by = "MULTI_ID", pt.size = 0.01, log = TRUE) + labs(title = sample_info$sample_id[i])
  print(HTO.vln.plots[[i]])
  dev.off()
  
  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("RidgePlot_nCountHTO_",sample_info$sample_id[i],".pdf")))
  HTO.ridge.plots[[i]] <- RidgePlot(seurat_pos[,!(seurat_pos@meta.data$MULTI_ID %in% c("Doublet"))], features = rownames(seurat_pos@assays$HTO@counts), group.by = "MULTI_ID", ncol = 2) + labs(title = sample_info$sample_id[i])
  print(HTO.ridge.plots[[i]])
  dev.off()
  
  # pdf(file = here::here("analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("FeatureScattert_",sample_info$sample_id[i],".pdf")))
  # HTO.scatter.plots[[i]] <- FeatureScatter(seurat_pos, feature1 = "hto_HTO-0251", feature2 = "hto_HTO-0257") + labs(title = sample_info$sample_id[i]) ##chose two that are in all samples (1 and 7)
  # print(HTO.scatter.plots[[i]])
  # dev.off()

  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("VlnPlot_nCount_RNA_",sample_info$sample_id[i],".pdf")))
  RNA.vln.plots[[i]] <- VlnPlot(seurat_pos, features = "nCount_RNA", pt.size = 0.1, log = TRUE) + labs(title = sample_info$sample_id[i])
  print(RNA.vln.plots[[i]])
  dev.off()
  
  message("Frequency of each hashtag identity including doublets and negatives\n",HTO.table <- plyr::count(seurat_pos@meta.data$MULTI_classification))
  write.csv(HTO.table, file = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_HTO_counts_preQC.csv")))
  
  #add doublet and negative counts to sample_info
  multi.id.tab <- plyr::count(seurat_pos@meta.data$MULTI_ID)

  # sample_info$doublets <- NA
  if(length(grep("Doublet",multi.id.tab$x)) > 0){
    sample_info$doublets[i] <- multi.id.tab[grep("Doublet",multi.id.tab$x),2]
  }
  sample_info$negatives[i] <- multi.id.tab[grep("Negative",multi.id.tab$x),2]

  # Subset only singlets. get rid of doublets and negatives
  seurat_pos <- seurat_pos %>% filter(!MULTI_ID %in% c("Negative","Doublet"))
  # seurat_pos <- subset(seurat_pos, idents %in% c("Negative", "Doublet"), invert = TRUE)
  sample_info$num.final.cells[i] <- dim(seurat_pos)[2]

  RNA.count.vlnplot[[i]] <- (VlnPlot(seurat_pos, features = "nCount_RNA",pt.size = 0) - VlnPlot(seurat_scran, features = "nCount_RNA",pt.size = 0)) / VlnPlot(seurat, features = "nCount_RNA",pt.size = 0,group.by = "orig.ident" )
  
  pdf(file = here::here("04_Analysis","plots","01_build_seurat",sample_info$sample_id[i],paste0("QCplots_after_DropCells",sample_info$sample_id[i],".pdf")))
  QCplots(seurat_pos)
  QCplots(seurat_scran)
  dev.off()

  QC.plots.after[[i]] <- QCplots(seurat_pos)
  #save both postive and neative seurat objects
  saveRDS(seurat_pos, file = here::here("04_Analysis","data_objects","01_build_seurat",sample_info$sample_id[i],paste0(sample_info$sample_id[i],"_pos.rds")))
  
  #merge with dropped cells dataframe
  removed.tab <- bind_rows(removed.tab, data.frame(Barcodes = rnascreenedDroppedCells, Reason = "01_build_seurat: RNA counts were either too high (doublet) or too low (empty droplet)"))
  removed.tab <- bind_rows(removed.tab, data.frame(Barcodes = dropcellsfxnDroppedCells, Reason = "01_build_seurat: Screened out by DropCells (dead or low quality cell- mitochondrial DNA and HK genes)"))
}

#save sample info table 
write.csv(sample_info, file = here::here("04_Analysis","data_objects","01_build_seurat","sample_info_stats.csv"))

#save grob lists 
# saveRDS(QC.plots.b4, file = here::here("analysis","plots","01_build_seurat","qc_b4_grobs_run2.rds")) 
# saveRDS(QC.plots.after, file = here::here("analysis","plots","01_build_seurat","qc_after_grobs_run2.rds"))
# saveRDS(background.plots, file = here::here("analysis","plots","01_build_seurat","background_grobs_run2.rds"))
# saveRDS(HTO.vln.plots,file = here::here("analysis","plots","01_build_seurat","HTO_vln_grobs_run2.rds"))
# saveRDS(HTO.ridge.plots, file = here::here("analysis","plots","01_build_seurat","HTO_ridge_grobs_run2.rds"))
# saveRDS(HTO.scatter.plots, file = here::here("analysis","plots","01_build_seurat","HTO_scatter_grobs_run2.rds"))
# saveRDS(RNA.vln.plots, file = here::here("analysis","plots","01_build_seurat","rna_vln_grobs_run2.rds"))

#plot grob lists
pdf(file = here::here("04_Analysis","plots","01_build_seurat","qc_b4_grobs.pdf"))
lapply(QC.plots.b4,plot)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","qc_after_grobs.pdf"))
lapply(QC.plots.after,plot)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","background_grobs.pdf"))
marrangeGrob(background.plots, ncol = 2, nrow = 2)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","HTO_vln_grobs.pdf"))
lapply(HTO.vln.plots, plot)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","HTO_ridge_grobs.pdf"))
lapply(HTO.ridge.plots,plot)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","HTO_scatter_grobs.pdf"))
marrangeGrob(HTO.scatter.plots, ncol = 2, nrow = 2)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","rna_vln_grobs.pdf"))
lapply(RNA.vln.plots, plot)
dev.off()

pdf(file = here::here("04_Analysis","plots","01_build_seurat","rna_count_vln_grobs.pdf"))
lapply(RNA.count.vlnplot, plot)
dev.off()

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()