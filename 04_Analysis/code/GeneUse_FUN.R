########### Plots for VDJ usage - completed ##################

GeneUse <- function(vdj, seurat.object = NULL, heavy.cutoff = 0, light.cutoff = 0, show.count = TRUE) { 
  
  H.names.index <- names(which(table(vdj$v_gene)/nrow(vdj) > heavy.cutoff))
  L.names.index <- names(which(table(vdj$v_gene.l)/nrow(vdj) > light.cutoff))
  
  
  if(show.count == TRUE) {
    counts <- geom_text(aes(label=..count..),stat="count",position=position_stack(0.5))
  } else {
    counts <- NULL
  }
  
  if (is.null(seurat.object)) { 
    p1 <- vdj %>% dplyr::select(v_gene, v_gene.l) %>%
      dplyr::filter(v_gene %in% H.names.index) %>%
      ggplot(aes(x = reorder(v_gene), y = ..count.. / nrow(vdj), fill = v_gene)) +
      geom_bar(position = "dodge2", alpha = 0.7, color = "black") + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.position = "none") + 
      ylab("Proportion") + 
      xlab("V Heavy Gene") + 
      counts + 
      scale_y_continuous(expand = c(0.0015, 0.0015))
    
    
    p2 <- vdj %>% dplyr::select(v_gene, v_gene.l) %>%
      dplyr::filter(v_gene.l %in% L.names.index) %>%
      ggplot(aes(x = v_gene.l, y = ..count.. / nrow(vdj), fill = v_gene.l)) +
      geom_bar(position = "dodge2", alpha = 0.7, color = "black") + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
            legend.position = "none") + 
      ylab("Proportion") + 
      xlab("V Light Gene") + 
      counts + 
      scale_y_continuous(expand = c(0.0015, 0.0015))
    
    grid.arrange(p1, p2, ncol = 2)
    
  } else {
    merger <- data.frame(barcode = colnames(seurat.object), HTO = as.factor(seurat.object@meta.data$HTO_classification))
    vdj$barcode <- gsub("-1|_1", "", vdj$barcode)
    vdj <- merge(vdj, merger, type = "left")
    grbs <- list()
    
    for ( i in 1:length(unique(vdj$HTO))) { 
      
      temp.HTO <- unique(vdj$HTO)[i]
      
      p1 <-  vdj %>% 
        dplyr::filter(v_gene %in% H.names.index & HTO == temp.HTO) %>%
        ggplot(aes(x = v_gene, y = ..count.. / sum(vdj$HTO == temp.HTO), fill = v_gene)) +
        geom_bar(position = "dodge2", alpha = 0.7, color = "black") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
              legend.position = "none") + 
        ylab("Proportion from HTO") + 
        xlab("V Heavy Gene") + 
        counts + ggtitle(paste(temp.HTO)) +
        scale_y_continuous(expand = c(0.0015, 0.0015))
      
      p2 <-   vdj %>%
        dplyr::filter(v_gene.l %in% L.names.index & HTO == temp.HTO) %>%
        ggplot(aes(x = v_gene.l, y = ..count.. / sum(vdj$HTO == temp.HTO), fill = v_gene.l)) +
        geom_bar(position = "dodge2", alpha = 0.7, color = "black") + 
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
              legend.position = "none") + 
        ylab("Proportion from HTO") + 
        xlab("V Light Gene") + 
        counts + ggtitle(paste(temp.HTO)) +
        scale_y_continuous(expand = c(0.0015, 0.0015))
      
      
      grbs[[i]] <- list(p1, p2)
      
    }
    
    grid.arrange(grobs = unlist(grbs, recursive = FALSE), ncol = 2)  
    
  }
}


###############################################