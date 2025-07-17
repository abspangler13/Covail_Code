######## Clonotyping - completed #######


Clonotype <- function(vdj, clonotype) { 
  
  
  clones2 <- clonotype %>% dplyr::select(SEQUENCE_ID, CLONE, V_CALL, D_CALL, J_CALL, C_CALL, SEQUENCE_INPUT, GERMLINE_IMGT, SEQUENCE_IMGT)
  clones2$SEQUENCE_ID <- gsub("_contig.*","", clones2$SEQUENCE_ID)
  
  file.names <- do.call(rbind, strsplit(as.character(clones2$SEQUENCE_ID),'_',fixed=TRUE))
  
  if(is.numeric(which(nchar(file.names[1,]) == 18))) {
    n <- which(nchar(file.names[1,]) == 18 )
  } else {
    n <- rev(order(nchar(file.names[1,])))[1] 
  }
  
  f.names <- file.names[,-n] 
  if (is.vector(f.names)) { f.names <- f.names
  } else {
    list.names <- list()
    for ( i in 1:ncol(f.names)) {
      
      list.names[[i]] <- f.names[,i]  
    }
    
    f.names <- gsub(" ", "_", do.call(paste, list.names))
  }
  
  clones3 <- data.frame(file = f.names, barcode = file.names[,n], clones2 %>% dplyr::select(-SEQUENCE_ID))
  
  
  vdj$barcode <- gsub("-1|_1", "", vdj$barcode)
  clones3$barcode <- gsub("-1|_1", "", clones3$barcode)
  final <- merge(vdj, clones3, by = "barcode", type = "left")
  msg <- table(final$file)
  message("barcodes found in :")
  message(paste0(capture.output(msg), collapse = "\n"))
  
  return(final)
  
}


#########################################