############### QC Plots - completed ###################


QCGenes <- function(files.list) {
  
  if (!is.list(files.list))
  { x <- list(files.list)
  } else {
    x <- files.list}
  
  house.keeper <- c("ACTB","B2M","HNRPLL","HPRT","PSMB2","PSMB4","PPIA","PRPS1","PRPS1L1",
                    "PRPS1L3","PRPS2","PRPSAP1","PRPSAP2","RPL10","RPL10A","RPL10L","RPL11",
                    "RPL12","RPL13","RPL14","RPL15","RPL17","RPL18","RPL19","RPL21","RPL22",
                    "RPL22L1","RPL23","RPL24","RPL26","RPL27","RPL28","RPL29","RPL3","RPL30",
                    "RPL32","RPL34","RPL35","RPL36","RPL37","RPL38","RPL39","RPL39L","RPL3L",
                    "RPL4","RPL41","RPL5","RPL6","RPL7","RPL7A","RPL7L1","RPL8","RPL9","RPLP0",
                    "RPLP1","RPLP2","RPS10","RPS11","RPS12","RPS13","RPS14","RPS15","RPS15A",
                    "RPS16","RPS17","RPS18","RPS19","RPS20","RPS21","RPS24","RPS25","RPS26",
                    "RPS27","RPS27A", "RPS27L","RPS28","RPS29","RPS3","RPS3A","RPS4X","RPS5",
                    "RPS6","RPS6KA1","RPS6KA2", "RPS6KA3","RPS6KA4","RPS6KA5","RPS6KA6","RPS6KB1",
                    "RPS6KB2","RPS6KC1","RPS6KL1","RPS7","RPS8","RPS9","RPSA","TRPS1","UBB")
  
  for (i in 1:length(x)) { 
    
    temp.counts <- x[[i]]@assays$RNA@counts
    
    mito.genes <- grep(pattern = "^MT-", x = rownames(x[[i]]), value = TRUE)
    
    percent.mito.temp <- Matrix::colSums(temp.counts[mito.genes, ])/Matrix::colSums(temp.counts)*100
    
    x[[i]] <- AddMetaData(x[[i]], metadata = percent.mito.temp, col.name = "percent.mito")
    
    hkgenes.found <- which(rownames(x[[i]]) %in% house.keeper)
    
    expressed.hkgenes <- Matrix::colSums(temp.counts[hkgenes.found,] > 0)
    
    
    x[[i]] <- AddMetaData(x[[i]], metadata = expressed.hkgenes, col.name = "n.exp.hkgenes")
    
  }
  
  avg <- cbind(
    unlist(
      lapply(x, function(x) { 
        mean(x@meta.data$percent.mito)
        }
        )
      ),
               
    unlist(
      lapply(x, function(x) { 
        mean(x@meta.data$n.exp.hkgene)
        }
        )
      ),
               
    unlist(
      lapply(x, function(x) { 
        mean(x@meta.data$nFeature_RNA)
        }
        )
      ),
               
    unlist(
      lapply(x, function(x) { 
        mean(x@meta.data$nCount_RNA)
        }
        )
      ), 
               
    unlist(
      lapply(x, function(x) { 
        nrow(x@meta.data)
        }
        )
      )
    )
  
  rownames(avg) <- unlist(lapply(x, function(x) { paste(x@project.name)}))
  colnames(avg) <- c("avg mito", "avg hk", "avg # genes", "avg counts", "# cells")
  message(paste0(capture.output(avg), collapse = "\n"))
  return(x)
  
}



########################################



# Warning messages:
#   1: In install.packages("tidyseurat") :
#   installation of package ‘SeuratObject’ had non-zero exit status
# 2: In install.packages("tidyseurat") :
#   installation of package ‘sctransform’ had non-zero exit status
# 3: In install.packages("tidyseurat") :
#   installation of package ‘Seurat’ had non-zero exit status
# 4: In install.packages("tidyseurat") :
#   installation of package ‘tidyseurat’ had non-zero exit status