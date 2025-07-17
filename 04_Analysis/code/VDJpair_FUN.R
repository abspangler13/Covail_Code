################# Pair function - completed ########################

VDJpair <- function(files, idents = NULL) {
  
  if (class(files) == "data.frame"){ 
    files <- list(files)
    } 
  else {
    files <- files }
  
  temp.file <- list()
  paired.files <- list()
  
  for (i in 1:length(files)) {
    
    temp.file <- files[[i]] %>% dplyr::filter(is_cell %in% c("TRUE","True") & full_length %in% c("TRUE","True")  & productive %in% c("TRUE","True")) %>% 
      dplyr::select(-is_cell, -contig_id, -high_confidence, -full_length, -productive, -raw_clonotype_id, -raw_consensus_id)
    
    temp.file.h <- temp.file %>% dplyr::filter(chain == "IGH")
    temp.file.k <- temp.file %>% dplyr::filter(chain == "IGK")
    temp.file.l <- temp.file %>% dplyr::filter(chain == "IGL")
    
    #remove duplicate barcodes
    temp.file.h <- temp.file.h[!(temp.file.h$barcode %in% (temp.file.h[duplicated(temp.file.h$barcode),]$barcode)),] 
    #refine this temp.file.h <- temp.file.h[!duplicated(temp.file.h$barcode),]
    temp.file.k <- temp.file.k[!(temp.file.k$barcode %in% (temp.file.k[duplicated(temp.file.k$barcode),]$barcode)),]
    temp.file.l <- temp.file.l[!(temp.file.l$barcode %in% (temp.file.l[duplicated(temp.file.l$barcode),]$barcode)),]
    
    
    temp.HK <- plyr::join(temp.file.h, temp.file.k, by = "barcode", type = "inner")
    temp.HL <- plyr::join(temp.file.h, temp.file.l, by = "barcode", type = "inner")
    
    colnames(temp.HK) <- c(colnames(temp.HK)[1:31], paste(colnames(temp.HK)[2:31], ".l", sep = ""))
    colnames(temp.HL) <- c(colnames(temp.HL)[1:31], paste(colnames(temp.HL)[2:31], ".l", sep = ""))
    
    temp.merge <- rbind(temp.HK, temp.HL)
    
    temp.merge <- temp.merge[!(temp.merge$barcode %in% temp.merge[duplicated(temp.merge$barcode),]$barcode),]
    
    if(length(idents) == length(files)) {
      temp.merge$orig <- paste(idents[i])
    } else {
      temp.merge$orig <- paste("file", i, sep = "_")
    }
    
    paired.files[[i]] <- temp.merge
  }
  
  rtrn <- (do.call(rbind, paired.files))
  return(rtrn %>% dplyr::select(-d_gene.l))
}
###############################################