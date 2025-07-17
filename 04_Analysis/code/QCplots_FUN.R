############### QC Plots - completed ###################


QCplots <- function(files.list) {
  
  if (!is.list(files.list))
  { x <- list(files.list)
  } else {
    x <- files.list}
  
  temp.grobs <- list()
  
  clrs <- c("wheat", "darkslateblue", "salmon", "green", "darkgray", "orange", "black")
  
  
  for (i in 1:length(x)) { 
    
    
    p1 <- x[[i]]@meta.data %>% 
      ggplot(aes(x = nFeature_RNA, y = nCount_RNA)) + 
      geom_point(color = clrs[i]) + 
      theme_classic() +
      theme(legend.position = "none") +
      ggtitle(paste(x[[i]]@project.name)) + 
      xlab("Number of Genes") + 
      ylab("Number of Transcripts")
    
    p2 <- x[[i]]@meta.data %>%
      ggplot(aes(x = n.exp.hkgenes, y = percent.mito)) + 
      geom_point(color = clrs[i]) + 
      theme_classic() +
      theme(legend.position = "none") +
      ggtitle(paste(x[[i]]@project.name)) + 
      xlab("Number of HK genes") +
      ylab("Percent mitochondrial")
    
    p3 <-  x[[i]]@meta.data %>%
      ggplot(aes(x = n.exp.hkgenes, color = "black")) + 
      geom_histogram(color = "black", bins = 35, fill = clrs[i]) + 
      theme_classic() +
      theme(legend.position = "none") +
      ggtitle(paste(x[[i]]@project.name)) + 
      xlab("Number of HK Genes") + 
      ylab("Count")
    
    temp.grobs[[i]] <- list(p1, p2, p3)
    
  }
  
  all.plots <- unlist(temp.grobs, recursive = FALSE)
  grid.arrange(grobs = all.plots, ncol = 3)
  
}



########################################
