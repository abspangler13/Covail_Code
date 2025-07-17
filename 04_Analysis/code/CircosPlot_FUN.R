################ Circos function - completed ############

CircosPlot <- function(vdj,  cutoff = 0) {
  
  x <- vdj
  
  for.circ <- unique(x[,c("v_gene",  "v_gene.l")])
  for.circ <- data.frame(for.circ, value = rep(0, nrow(for.circ)))
  
  
  genes <- x %>%dplyr:: select(v_gene,  v_gene.l)
  
  for (i in 1:nrow(genes)) { 
    
    temp.H <- genes[i,1]
    temp.L <- genes[i,2]
    
    a <- which(for.circ$v_gene == temp.H)
    b <- which(for.circ$v_gene.l == temp.L)
    c <- intersect(a, b)
    
    for.circ[c,3] <- for.circ[c,3] + 1
    
  }
  
  for.circ
  c <- for.circ[for.circ$value > cutoff,]
  
  par(mfrow = c(1,1))
  par(mar = c(4,4,4,4))
  
  circos.clear()
  circos.par(cell.padding=c(1,1,1,1), track.margin=c(0.0,0.15), start.degree = 0, gap.degree = 1, canvas.ylim = c(-1.5,1.5))
  chordDiagram(c, annotationTrack = "grid")
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(
      x = CELL_META$xcenter, 
      y = CELL_META$ylim[1] + 1, 
      labels = CELL_META$sector.index, 
      facing = "clockwise", 
      niceFacing = TRUE,
      adj = c(-0.2, 0.1),
      cex = 0.9)}, bg.border = NA)
  
  
} 
#################################################
