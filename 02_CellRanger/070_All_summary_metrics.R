library(ggplot2)
library(patchwork)
library(tidyverse)

files <- list.files("CellRanger_run5",pattern = "(-VDJ|-CSO|-GEX)$")
vdj.dat <- data.frame(matrix(ncol = 33, nrow = 0))
cso.dat <- data.frame(matrix(ncol = 15, nrow = 0))
gex.dat <- data.frame(matrix(ncol = 19, nrow = 0))

for(i in 1:length(files)){
  x <- read.csv(here::here("CellRanger_run5",files[i],"outs","metrics_summary.csv"))
  if(ncol(x)==15){
    cso.dat <- rbind(cso.dat,x)
  }
  if(ncol(x)==19){
    gex.dat <- rbind(gex.dat,x)
  }
  if(ncol(x)==33){
    vdj.dat <- rbind(vdj.dat,x)
  }
}

rownames(vdj.dat) <- files[grep("-VDJ",files)]
rownames(cso.dat) <- files[grep("-CSO",files)]
rownames(gex.dat) <- files[grep("-GEX",files)]

vdj.dat$Estimated.Number.of.Cells <- as.numeric(gsub(",","",vdj.dat$Estimated.Number.of.Cells))
vdj.dat$Mean.Read.Pairs.per.Cell  <- as.numeric(gsub(",","",vdj.dat$Mean.Read.Pairs.per.Cell ))
vdj.dat$Number.of.Cells.With.Productive.V.J.Spanning.Pair <- as.numeric(gsub(",","",vdj.dat$Number.of.Cells.With.Productive.V.J.Spanning.Pair))
vdj.dat$Number.of.Read.Pairs <- as.numeric(gsub(",","",vdj.dat$Number.of.Read.Pairs))

cso.dat$Estimated.Number.of.Cells <- as.numeric(gsub(",","",cso.dat$Estimated.Number.of.Cells))
cso.dat$Antibody..Mean.Reads.per.Cell  <- as.numeric(gsub(",","",cso.dat$Antibody..Mean.Reads.per.Cell))

gex.dat$Estimated.Number.of.Cells <- as.numeric(gsub(",","",gex.dat$Estimated.Number.of.Cells))
gex.dat$Mean.Reads.per.Cell  <- as.numeric(gsub(",","",gex.dat$Mean.Reads.per.Cell))


dat.obs <- list(vdj.dat,cso.dat,gex.dat)
names(dat.obs) <- c("vdj","cso","gex")
for(i in 1:length(dat.obs)){
  dat.obs[[i]]$library <- rownames(dat.obs[[i]])
  dat.obs[[i]]$cell_type[grep("PB",dat.obs[[i]]$library)] <- "PB"
  dat.obs[[i]]$cell_type[grep("M",dat.obs[[i]]$library)] <- "Mem"
  dat.obs[[i]]$sample <- gsub('.{4}$', '', dat.obs[[i]]$library)
  dat.obs[[i]]$library_type <- str_sub(dat.obs[[i]]$library,-3,-1)
  dat.obs[[i]]$donor <- str_sub(dat.obs[[i]]$library,7,9)
  dat.obs[[i]] <- dat.obs[[i]] %>% arrange(cell_type,donor)
  write.csv(dat.obs[[i]],file = here::here("CellRanger_run5",paste0(names(dat.obs[i]),"_summary_metrics.csv")))
}

###make plots
v1<-ggplot(data=dat.obs$vdj, aes(x=fct_reorder(sample,cell_type), y=Estimated.Number.of.Cells, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("VDJ Mean Number of Cells")

v2<-ggplot(data=dat.obs$vdj, aes(x=fct_reorder(sample,cell_type), y=Mean.Read.Pairs.per.Cell, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("VDJ Mean Reads per Cell")

c1<-ggplot(data=dat.obs$cso, aes(x=fct_reorder(sample,cell_type), y=Estimated.Number.of.Cells, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("CSO Mean Number of Cells")

c2<-ggplot(data=dat.obs$cso, aes(x=fct_reorder(sample,cell_type), y=Antibody..Mean.Reads.per.Cell, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("CSO Mean Reads per Cell")

g1<-ggplot(data=dat.obs$gex, aes(x=fct_reorder(sample,cell_type), y=Estimated.Number.of.Cells, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("GEX Mean Number of Cells")

g2<-ggplot(data=dat.obs$gex, aes(x=fct_reorder(sample,cell_type), y=Mean.Reads.per.Cell, fill = cell_type)) +
  geom_bar(stat="identity") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("GEX Mean Reads per Cell")

pdf(file = here::here("CellRanger_run5","QC_metrics_plot.pdf"), width = 16, height = 16)
(v1 + v2) / (c1 + c2) / (g1 + g2) 
dev.off()


v <- dat.obs$vdj %>% select(Estimated.Number.of.Cells,Mean.Read.Pairs.per.Cell,library,cell_type,library_type,sample,donor) %>% rename(Mean.Reads.per.Cell = Mean.Read.Pairs.per.Cell)
c <- dat.obs$cso %>% select(Estimated.Number.of.Cells,Antibody..Mean.Reads.per.Cell,library,cell_type,library_type,sample,donor) %>% rename(Mean.Reads.per.Cell = Antibody..Mean.Reads.per.Cell)
g <- dat.obs$gex %>% select(Estimated.Number.of.Cells,Mean.Reads.per.Cell,library,cell_type,library_type,sample,donor) 

All.dat <- bind_rows(v,c,g)

pdf(file = here::here("CellRanger_run5","Cell_efficiency.pdf"))
ggplot(data=All.dat, aes(x=fct_reorder(sample, cell_type), y=Estimated.Number.of.Cells, fill=library_type)) +
geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggtitle("Cell Number Efficiency")
dev.off()