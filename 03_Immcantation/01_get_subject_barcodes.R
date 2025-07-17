#Adapted from Abby's code
### Get sample barcodes for MemA and MemB pools
library(Seurat)
library(tidyseurat)
library(phylotools)
library(here)

#load in the seurat object containing GEX and CSO data so that we can pull out relevant barcodes
seuObj <- readRDS(file = here::here("04_Analysis","data_objects","01_build_seurat","MergedSeuratObject_p.rds"))

meta.data <- seuObj@meta.data %>% select(Subject)
meta.data$barcode <- rownames(meta.data)

#we will need to be careful in how we approach this
#First: pull gel bead barcodes per pool that correspond to subject so that we can demultiplex between two people per pool

#Second: Write a csv of these for safe keeping later

#Third: load in the fasta files, pull the sequences belonging to individual donors, and save them as separate FASTAs

#Fourth: if possible, I think we should append the pool names at this point. I don't know if it really matters adding the subject name simultaneously.
#When we add Immcantation data to the seurat object, we only really need the pool appended to contig name. This would eliminate the need to do the 00 script

#####
## "COV_09_13_2023_P1
s5351564848 <- meta.data %>% filter(Subject == 5351564848) %>% select(barcode)

write.csv(s5351564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5351564848_barcodes.csv"))

COV_09_13_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-13-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_13_2023_P1$seq.name <- paste0("COV_09_13_2023_P1_", COV_09_13_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5351564848
indices <- unlist(sapply(s5351564848$barcode, function(x){
  grep(x,COV_09_13_2023_P1$seq.name)
}))
fasta_5351564848 <- COV_09_13_2023_P1[indices,]
dat2fasta(fasta_5351564848, outfile = here::here("03_Immcantation","s5351564848_filtered_contig.fasta"))

rm(COV_09_13_2023_P1) #so I don't get confused

#####

#####
## "COV_09_19_2023_P1"
s5457484948 <- meta.data %>% filter(Subject == 5457484948) %>% select(barcode)
s5048574848 <- meta.data %>% filter(Subject == 5048574848) %>% select(barcode)

write.csv(s5457484948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5457484948_barcodes.csv"))
write.csv(s5048574848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5048574848_barcodes.csv"))

COV_09_19_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-19-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_19_2023_P1$seq.name <- paste0("COV_09_19_2023_P1_", COV_09_19_2023_P1$seq.name)

#5457484948
indices <- unlist(sapply(s5457484948$barcode, function(x){
  grep(x,COV_09_19_2023_P1$seq.name)
}))
fasta_5457484948 <- COV_09_19_2023_P1[indices,]
dat2fasta(fasta_5457484948, outfile = here::here("03_Immcantation","s5457484948_filtered_contig.fasta"))

#5048574848
indices <- unlist(sapply(s5048574848$barcode, function(x){
  grep(x,COV_09_19_2023_P1$seq.name)
}))
fasta_s5048574848 <- COV_09_19_2023_P1[indices,]
dat2fasta(fasta_s5048574848, outfile = here::here("03_Immcantation","s5048574848_filtered_contig.fasta"))

rm(COV_09_19_2023_P1)
#####

#####
# "COV_09_19_2023_P2/P3"
s4951574848 <- meta.data %>% filter(Subject == 4951574848) %>% select(barcode)
s5357484948 <- meta.data %>% filter(Subject == 5357484948) %>% select(barcode)
write.csv(s4951574848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4951574848_barcodes.csv"))
write.csv(s5357484948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5357484948_barcodes.csv"))

COV_09_19_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-19-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_19_2023_P2$seq.name <- paste0("COV_09_19_2023_P2_", COV_09_19_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_19_2023_P3 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-19-2023-VDJ-Pool3","outs","filtered_contig.fasta"))
COV_09_19_2023_P3$seq.name <- paste0("COV_09_19_2023_P3_", COV_09_19_2023_P3$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_19_2023P2P3 <- rbind(COV_09_19_2023_P2, COV_09_19_2023_P3)

#4951574848
indices <- unlist(sapply(s4951574848$barcode, function(x){
  grep(x,COV_09_19_2023P2P3$seq.name)
}))
fasta_4951574848 <- COV_09_19_2023P2P3[indices,]
dat2fasta(fasta_4951574848, outfile = here::here("03_Immcantation","s4951574848_filtered_contig.fasta"))

#5357484948
indices <- unlist(sapply(s5357484948$barcode, function(x){
  grep(x,COV_09_19_2023P2P3$seq.name)
}))
fasta_5357484948 <- COV_09_19_2023P2P3[indices,]
dat2fasta(fasta_5357484948, outfile = here::here("03_Immcantation","s5357484948_filtered_contig.fasta"))

rm(COV_09_19_2023_P2)
rm(COV_09_19_2023_P3)
rm(COV_09_19_2023P2P3)
#####

#####
# "COV_09_20_2023_P1"
s4954554848 <- meta.data %>% filter(Subject == 4954554848) %>% select(barcode)
s4957484948 <- meta.data %>% filter(Subject == 4957484948) %>% select(barcode)
write.csv(s4954554848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4954554848_barcodes.csv"))
write.csv(s4957484948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4957484948_barcodes.csv"))

COV_09_20_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_20_2023_P1$seq.name <- paste0("COV_09_20_2023_P1_", COV_09_20_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#4954554848
indices <- unlist(sapply(s4954554848$barcode, function(x){
  grep(x,COV_09_20_2023_P1$seq.name)
}))
fasta_4954554848 <- COV_09_20_2023_P1[indices,]
dat2fasta(fasta_4954554848, outfile = here::here("03_Immcantation","s4954554848_filtered_contig.fasta"))

#4957484948
indices <- unlist(sapply(s4957484948$barcode, function(x){
  grep(x,COV_09_20_2023_P1$seq.name)
}))
fasta_4957484948 <- COV_09_20_2023_P1[indices,]
dat2fasta(fasta_4957484948, outfile = here::here("03_Immcantation","s4957484948_filtered_contig.fasta"))

rm(COV_09_20_2023_P1)
#####

#####
# "COV_09_20_2023_P2/P3"
s4848544848 <- meta.data %>% filter(Subject == 4848544848) %>% select(barcode)
s5053564848 <- meta.data %>% filter(Subject == 5053564848) %>% select(barcode)
write.csv(s4848544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4848544848_barcodes.csv"))
write.csv(s5053564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5053564848_barcodes.csv"))

COV_09_20_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_20_2023_P2$seq.name <- paste0("COV_09_20_2023_P2_", COV_09_20_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_20_2023_P3 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool3","outs","filtered_contig.fasta"))
COV_09_20_2023_P3$seq.name <- paste0("COV_09_20_2023_P3_", COV_09_20_2023_P3$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_20_2023P2P3 <- rbind(COV_09_20_2023_P2, COV_09_20_2023_P3)

#4848544848
indices <- unlist(sapply(s4848544848$barcode, function(x){
  grep(x,COV_09_20_2023P2P3$seq.name)
}))
fasta_4848544848 <- COV_09_20_2023P2P3[indices,]
dat2fasta(fasta_4848544848, outfile = here::here("03_Immcantation","s4848544848_filtered_contig.fasta"))

#5053564848 
indices <- unlist(sapply(s5053564848 $barcode, function(x){
  grep(x,COV_09_20_2023P2P3$seq.name)
}))
fasta_5053564848  <- COV_09_20_2023P2P3[indices,]
dat2fasta(fasta_5053564848, outfile = here::here("03_Immcantation","s5053564848_filtered_contig.fasta"))

rm(COV_09_20_2023_P2)
rm(COV_09_20_2023_P3)
rm(COV_09_20_2023P2P3)
#####

#####
# "COV_09_21_2023_P1"
s4955534848 <- meta.data %>% filter(Subject == 4955534848) %>% select(barcode)
s4957564848 <- meta.data %>% filter(Subject == 4957564848) %>% select(barcode)
write.csv(s4955534848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4955534848_barcodes.csv"))
write.csv(s4957564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4957564848_barcodes.csv"))

COV_09_21_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-21-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_21_2023_P1$seq.name <- paste0("COV_09_21_2023_P1_", COV_09_21_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#4955534848
indices <- unlist(sapply(s4955534848$barcode, function(x){
  grep(x,COV_09_21_2023_P1$seq.name)
}))
fasta_4955534848 <- COV_09_21_2023_P1[indices,]
dat2fasta(fasta_4955534848, outfile = here::here("03_Immcantation","s4955534848_filtered_contig.fasta"))

#s4957564848
indices <- unlist(sapply(s4957564848$barcode, function(x){
  grep(x,COV_09_21_2023_P1$seq.name)
}))
fasta_s4957564848 <- COV_09_21_2023_P1[indices,]
dat2fasta(fasta_s4957564848, outfile = here::here("03_Immcantation","s4957564848_filtered_contig.fasta"))

rm(COV_09_21_2023_P1) #so I don't get confused
#####

#####
# "COV_09_21_2023_P2"
s4953494948 <- meta.data %>% filter(Subject == 4953494948) %>% select(barcode)
s5553564848 <- meta.data %>% filter(Subject == 5553564848) %>% select(barcode)
write.csv(s4953494948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4953494948_barcodes.csv"))
write.csv(s5553564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5553564848_barcodes.csv"))

COV_09_21_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-21-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_21_2023_P2$seq.name <- paste0("COV_09_21_2023_P2_", COV_09_21_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#4953494948
indices <- unlist(sapply(s4953494948$barcode, function(x){
  grep(x,COV_09_21_2023_P2$seq.name)
}))
fasta_4953494948 <- COV_09_21_2023_P2[indices,]
dat2fasta(fasta_4953494948, outfile = here::here("03_Immcantation","s4953494948_filtered_contig.fasta"))

#5553564848
indices <- unlist(sapply(s5553564848$barcode, function(x){
  grep(x,COV_09_21_2023_P2$seq.name)
}))
fasta_s5553564848 <- COV_09_21_2023_P2[indices,]
dat2fasta(fasta_s5553564848, outfile = here::here("03_Immcantation","s5553564848_filtered_contig.fasta"))

rm(COV_09_21_2023_P2)
#####

#####
# "COV_09_26_2023_P1"
s5150564848 <- meta.data %>% filter(Subject == 5150564848) %>% select(barcode)
s5557544848 <- meta.data %>% filter(Subject == 5557544848) %>% select(barcode)
write.csv(s5150564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5150564848_barcodes.csv"))
write.csv(s5557544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5557544848_barcodes.csv"))

COV_09_26_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-26-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_26_2023_P1$seq.name <- paste0("COV_09_26_2023_P1_", COV_09_26_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5150564848
indices <- unlist(sapply(s5150564848$barcode, function(x){
  grep(x,COV_09_26_2023_P1$seq.name)
}))
fasta_5150564848 <- COV_09_26_2023_P1[indices,]
dat2fasta(fasta_5150564848, outfile = here::here("03_Immcantation","s5150564848_filtered_contig.fasta"))

#5557544848
indices <- unlist(sapply(s5557544848$barcode, function(x){
  grep(x,COV_09_26_2023_P1$seq.name)
}))
fasta_s5557544848 <- COV_09_26_2023_P1[indices,]
dat2fasta(fasta_s5557544848, outfile = here::here("03_Immcantation","s5557544848_filtered_contig.fasta"))

rm(COV_09_26_2023_P1)
#####

#####
# "COV_09_26_2023_P2"
s5050484948 <- meta.data %>% filter(Subject == 5050484948) %>% select(barcode)
s5048544848 <- meta.data %>% filter(Subject == 5048544848) %>% select(barcode)
write.csv(s5050484948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5050484948_barcodes.csv"))
write.csv(s5048544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5048544848_barcodes.csv"))

COV_09_26_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-26-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_26_2023_P2$seq.name <- paste0("COV_09_26_2023_P2_", COV_09_26_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5050484948
indices <- unlist(sapply(s5050484948$barcode, function(x){
  grep(x,COV_09_26_2023_P2$seq.name)
}))
fasta_5050484948 <- COV_09_26_2023_P2[indices,]
dat2fasta(fasta_5050484948, outfile = here::here("03_Immcantation","s5050484948_filtered_contig.fasta"))

#5048544848
indices <- unlist(sapply(s5048544848$barcode, function(x){
  grep(x,COV_09_26_2023_P2$seq.name)
}))
fasta_s5048544848 <- COV_09_26_2023_P2[indices,]
dat2fasta(fasta_s5048544848, outfile = here::here("03_Immcantation","s5048544848_filtered_contig.fasta"))

rm(COV_09_26_2023_P2)
#####

#####
# "COV_09_27_2023_P1"
s5750564848 <- meta.data %>% filter(Subject == 5750564848) %>% select(barcode)
s5755544848 <- meta.data %>% filter(Subject == 5755544848) %>% select(barcode)
write.csv(s5750564848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5750564848_barcodes.csv"))
write.csv(s5755544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5755544848_barcodes.csv"))

COV_09_27_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-27-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_27_2023_P1$seq.name <- paste0("COV_09_27_2023_P1_", COV_09_27_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5750564848
indices <- unlist(sapply(s5750564848$barcode, function(x){
  grep(x,COV_09_27_2023_P1$seq.name)
}))
fasta_5750564848 <- COV_09_27_2023_P1[indices,]
dat2fasta(fasta_5750564848, outfile = here::here("03_Immcantation","s5750564848_filtered_contig.fasta"))

#5755544848
indices <- unlist(sapply(s5755544848$barcode, function(x){
  grep(x,COV_09_27_2023_P1$seq.name)
}))
fasta_s5755544848 <- COV_09_27_2023_P1[indices,]
dat2fasta(fasta_s5755544848, outfile = here::here("03_Immcantation","s5755544848_filtered_contig.fasta"))

rm(COV_09_27_2023_P1)
#####

#####
# "COV_09_27_2023_P2"
s4950544848 <- meta.data %>% filter(Subject == 4950544848) %>% select(barcode)
s5249544848 <- meta.data %>% filter(Subject == 5249544848) %>% select(barcode)
write.csv(s4950544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4950544848_barcodes.csv"))
write.csv(s5249544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5249544848_barcodes.csv"))

COV_09_27_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-27-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_27_2023_P2$seq.name <- paste0("COV_09_27_2023_P2_", COV_09_27_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#4950544848
indices <- unlist(sapply(s4950544848$barcode, function(x){
  grep(x,COV_09_27_2023_P2$seq.name)
}))
fasta_4950544848 <- COV_09_27_2023_P2[indices,]
dat2fasta(fasta_4950544848, outfile = here::here("03_Immcantation","s4950544848_filtered_contig.fasta"))

#5249544848
indices <- unlist(sapply(s5249544848$barcode, function(x){
  grep(x,COV_09_27_2023_P2$seq.name)
}))
fasta_s5249544848 <- COV_09_27_2023_P2[indices,]
dat2fasta(fasta_s5249544848, outfile = here::here("03_Immcantation","s5249544848_filtered_contig.fasta"))

rm(COV_09_27_2023_P2)
#####

#####
# "COV_09_28_2023_P1/P3"
s5054574848 <- meta.data %>% filter(Subject == 5054574848) %>% select(barcode)
s5456544848 <- meta.data %>% filter(Subject == 5456544848) %>% select(barcode)
write.csv(s5054574848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5054574848_barcodes.csv"))
write.csv(s5456544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5456544848_barcodes.csv"))

COV_09_28_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_09_28_2023_P1$seq.name <- paste0("COV_09_28_2023_P1_", COV_09_28_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_28_2023_P3 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool3","outs","filtered_contig.fasta"))
COV_09_28_2023_P3$seq.name <- paste0("COV_09_28_2023_P3_", COV_09_28_2023_P3$seq.name) #go ahead and paste sequence name so we don't need the 01 script

COV_09_28_2023P1P3 <- rbind(COV_09_28_2023_P1, COV_09_28_2023_P3)

#5054574848
indices <- unlist(sapply(s5054574848$barcode, function(x){
  grep(x,COV_09_28_2023P1P3$seq.name)
}))
fasta_5054574848 <- COV_09_28_2023P1P3[indices,]
dat2fasta(fasta_5054574848, outfile = here::here("03_Immcantation","s5054574848_filtered_contig.fasta"))

#5456544848
indices <- unlist(sapply(s5456544848$barcode, function(x){
  grep(x,COV_09_28_2023P1P3$seq.name)
}))
fasta_5456544848 <- COV_09_28_2023P1P3[indices,]
dat2fasta(fasta_5456544848, outfile = here::here("03_Immcantation","s5456544848_filtered_contig.fasta"))

rm(COV_09_28_2023_P1)
rm(COV_09_28_2023_P3)
rm(COV_09_28_2023P2P3)
#####

#####
# "COV_09_28_2023_P2"
s5356534848 <- meta.data %>% filter(Subject == 5356534848) %>% select(barcode)
s5556484948 <- meta.data %>% filter(Subject == 5556484948) %>% select(barcode)
write.csv(s5356534848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5356534848_barcodes.csv"))
write.csv(s5556484948, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5556484948_barcodes.csv"))

COV_09_28_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_09_28_2023_P2$seq.name <- paste0("COV_09_28_2023_P2_", COV_09_28_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5356534848
indices <- unlist(sapply(s5356534848$barcode, function(x){
  grep(x,COV_09_28_2023_P2$seq.name)
}))
fasta_5356534848 <- COV_09_28_2023_P2[indices,]
dat2fasta(fasta_5356534848, outfile = here::here("03_Immcantation","s5356534848_filtered_contig.fasta"))

#5556484948
indices <- unlist(sapply(s5556484948$barcode, function(x){
  grep(x,COV_09_28_2023_P2$seq.name)
}))
fasta_s5556484948 <- COV_09_28_2023_P2[indices,]
dat2fasta(fasta_s5556484948, outfile = here::here("03_Immcantation","s5556484948_filtered_contig.fasta"))

rm(COV_09_28_2023_P2)
#####

#####
# "COV_10_05_2023_P1"
s5553554848 <- meta.data %>% filter(Subject == 5553554848) %>% select(barcode)
write.csv(s5553554848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5553554848_barcodes.csv"))

COV_10_05_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-10-05-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_10_05_2023_P1$seq.name <- paste0("COV_10_05_2023_P1_", COV_10_05_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5553554848
indices <- unlist(sapply(s5553554848$barcode, function(x){
  grep(x,COV_10_05_2023_P1$seq.name)
}))
fasta_5553554848 <- COV_10_05_2023_P1[indices,]
dat2fasta(fasta_5553554848, outfile = here::here("03_Immcantation","s5553554848_filtered_contig.fasta"))

rm(COV_10_05_2023_P1)
#####

#####
# "COV_10_06_2023_P1"
s5657574848 <- meta.data %>% filter(Subject == 5657574848) %>% select(barcode)
write.csv(s5657574848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5657574848_barcodes.csv"))

COV_10_06_2023_P1 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-10-06-2023-VDJ-Pool1","outs","filtered_contig.fasta"))
COV_10_06_2023_P1$seq.name <- paste0("COV_10_06_2023_P1_", COV_10_06_2023_P1$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#5657574848
indices <- unlist(sapply(s5657574848$barcode, function(x){
  grep(x,COV_10_06_2023_P1$seq.name)
}))
fasta_5657574848 <- COV_10_06_2023_P1[indices,]
dat2fasta(fasta_5657574848, outfile = here::here("03_Immcantation","s5657574848_filtered_contig.fasta"))

rm(COV_10_06_2023_P1)
#####

#####
# "COV_10_06_2023_P2"
s4856554848 <- meta.data %>% filter(Subject == 4856554848) %>% select(barcode)
s5653544848 <- meta.data %>% filter(Subject == 5653544848) %>% select(barcode)
write.csv(s4856554848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s4856554848_barcodes.csv"))
write.csv(s5653544848, file = here::here("03_Immcantation","CorrectedFASTAFiles","s5653544848_barcodes.csv"))

COV_10_06_2023_P2 <- phylotools::read.fasta(file = here::here("02_CellRanger","VRC-COVAIL-10-06-2023-VDJ-Pool2","outs","filtered_contig.fasta"))
COV_10_06_2023_P2$seq.name <- paste0("COV_10_06_2023_P2_", COV_10_06_2023_P2$seq.name) #go ahead and paste sequence name so we don't need the 01 script

#4856554848
indices <- unlist(sapply(s4856554848$barcode, function(x){
  grep(x,COV_10_06_2023_P2$seq.name)
}))
fasta_4856554848 <- COV_10_06_2023_P2[indices,]
dat2fasta(fasta_4856554848, outfile = here::here("03_Immcantation","s4856554848_filtered_contig.fasta"))

#5653544848
indices <- unlist(sapply(s5653544848$barcode, function(x){
  grep(x,COV_10_06_2023_P2$seq.name)
}))
fasta_s5653544848 <- COV_10_06_2023_P2[indices,]
dat2fasta(fasta_s5653544848, outfile = here::here("03_Immcantation","s5653544848_filtered_contig.fasta"))

rm(COV_10_06_2023_P2)
#####
