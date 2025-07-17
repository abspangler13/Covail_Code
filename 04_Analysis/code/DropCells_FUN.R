################## Drop Cells - conpleted #############


DropCells <- function(data, min.features, max.mito, min.hk, max.count, min.count, path.out = 0) {
  
  x <- data
  
  index <- data@meta.data %>% dplyr::mutate(barcodes = rownames(x@meta.data)) %>% dplyr::filter(nCount_RNA >= min.count &
                                                                                                  nCount_RNA <= max.count & max.mito >= percent.mito & n.exp.hkgenes >= min.hk)
  barcodes.pass <- index$barcodes
  freq <- ((ncol(x) - length(barcodes.pass))/ncol(x))*100
  message("Dropped ", paste(ncol(x) - length(barcodes.pass), "Cells (", format(round(freq, 4), nsmall = 4),"% )"))
  
  
  barcodes.fail <- rownames(x@meta.data)[!(rownames(x@meta.data) %in% barcodes.pass)]
  
  fail.table <- x@meta.data[barcodes.fail,]
  
  reason.failed <- data.frame(0, 0, 0, 0, 0)
  colnames(reason.failed) <- c(paste("RNA >", max.count), paste("RNA <", min.count), paste("Genes <", min.features), paste("Mito >", max.mito, "%"), paste("hk genes <", min.hk))
  
  for (i in 1:nrow(fail.table)) { 
    temp.row <- fail.table[i,]
    col1 <- temp.row$nCount_RNA < min.count
    col2 <- temp.row$nCount_RNA > max.count
    col3 <- temp.row$nFeature_RNA < min.features
    col4 <- temp.row$percent.mito > max.mito
    col5 <- temp.row$n.exp.hkgenes < min.hk
    
    reason.failed[i,1:5] <- c(col1, col2, col3, col4, col5)
  }
  
  outs <- cbind(fail.table, reason.failed)
  
  means <- vector()
  
  for (i in 1:ncol(outs)) { 
    if(class(outs[,i]) %in% c("numeric", "integer")) {
      means[i] <- mean(outs[,i])
    } else {
      means[i] <- NA 
    }
  }
  
  means <- suppressWarnings(as.numeric(format(round(means, 4), nsmall = 4)))
  outs2 <- rbind(outs, AVG = means)
  if (is.null(path.out)) { write.csv(x = outs2, paste(getwd(),"/",paste("QC_drop_", x@project.name,".csv", sep = ""), sep =""))}
  else if (path.out == 0) { 
    message("No QC file generated")
  } else {
    # write.csv(x = outs2, paste(fp,"/",paste("QC_drop_", x@project.name,".csv", sep = ""), sep =""))
    write.csv(x = outs2, file = path.out)
  }
  return(subset(x, cells = barcodes.pass))
}

#############################################
