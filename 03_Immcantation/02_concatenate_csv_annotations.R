#Making proper csv annotation files

#In order to do this, we'll need to grab the barcodes we want this time- we wrote this from the previous step
#After we load the barcodes, we'll load the filtered contig annotations and grab ones that correspond to barcodes we want


#####
#load the barcode files
s5351564848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5351564848_barcodes.csv"))
s5457484948 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5457484948_barcodes.csv"))
s5048574848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5048574848_barcodes.csv"))
s4951574848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4951574848_barcodes.csv"))
s5357484948 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5357484948_barcodes.csv"))
s4954554848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s4954554848_barcodes.csv"))
s4957484948 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4957484948_barcodes.csv"))
s4848544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4848544848_barcodes.csv"))
s5053564848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5053564848_barcodes.csv"))
s4955534848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4955534848_barcodes.csv"))
s4957564848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4957564848_barcodes.csv"))
s4953494948 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4953494948_barcodes.csv"))
s5553564848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5553564848_barcodes.csv"))
s5150564848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5150564848_barcodes.csv"))
s5557544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5557544848_barcodes.csv"))
s5050484948 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5050484948_barcodes.csv"))
s5048544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5048544848_barcodes.csv"))
s5750564848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5750564848_barcodes.csv"))
s5755544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5755544848_barcodes.csv"))
s4950544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4950544848_barcodes.csv"))
s5249544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5249544848_barcodes.csv"))
s5054574848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5054574848_barcodes.csv"))
s5456544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5456544848_barcodes.csv"))
s5356534848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5356534848_barcodes.csv"))
s5556484948 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5556484948_barcodes.csv"))
s5553554848 <- read.csv(file = here::here("03_Immcantation","CorrectedFASTAFiles","s5553554848_barcodes.csv"))
s5657574848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5657574848_barcodes.csv"))
s4856554848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s4856554848_barcodes.csv"))
s5653544848 <- read.csv( file = here::here("03_Immcantation","CorrectedFASTAFiles","s5653544848_barcodes.csv"))
#####

#####
#Begin separating .csv by actual subject
#####
##"COV_09_13_2023_P1"
COV_09_13_2023_P1 <- read.csv(here::here("02_CellRanger", "VRC-COVAIL-09-13-2023-VDJ-Pool1", "outs", "filtered_contig_annotations.csv"))
COV_09_13_2023_P1$barcode <- paste0("COV_09_13_2023_P1_", COV_09_13_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_13_2023_P1$contig_id <- paste0("COV_09_13_2023_P1_", COV_09_13_2023_P1$contig_id)

#s5351564848
indices <- unlist(sapply(s5351564848$barcode, function(x){
  grep(x,COV_09_13_2023_P1$barcode)
}))
csv_s5351564848 <- COV_09_13_2023_P1[indices,]
write.csv(csv_s5351564848, file=here::here("03_Immcantation","s5351564848_filtered_contig_annotations.csv"))

rm(COV_09_13_2023_P1)
####

## "COV_09_19_2023_P1"
COV_09_19_2023_P1 <- read.csv(here::here("02_CellRanger", "VRC-COVAIL-09-19-2023-VDJ-Pool1", "outs", "filtered_contig_annotations.csv"))
COV_09_19_2023_P1$barcode <- paste0("COV_09_19_2023_P1_", COV_09_19_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_19_2023_P1$contig_id <- paste0("COV_09_19_2023_P1_", COV_09_19_2023_P1$contig_id)

#5457484948
indices <- unlist(sapply(s5457484948$barcode, function(x){
  grep(x,COV_09_19_2023_P1$barcode)
}))
csv_5457484948 <- COV_09_19_2023_P1[indices,]
write.csv(csv_5457484948, file=here::here("03_Immcantation","s5457484948_filtered_contig_annotations.csv"))

#5048574848
indices <- unlist(sapply(s5048574848$barcode, function(x){
  grep(x,COV_09_19_2023_P1$barcode)
}))
csv_s5048574848 <- COV_09_19_2023_P1[indices,]
write.csv(csv_s5048574848, file = here::here("03_Immcantation","s5048574848_filtered_contig_annotations.csv"))

rm(COV_09_19_2023_P1) #so I don't get confused
#####

#####
# "COV_09_19_2023_P2/P3"
COV_09_19_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-19-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_19_2023_P2$barcode <- paste0("COV_09_19_2023_P2_", COV_09_19_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_19_2023_P2$contig_id <- paste0("COV_09_19_2023_P2_", COV_09_19_2023_P2$contig_id)

COV_09_19_2023_P3 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-19-2023-VDJ-Pool3","outs","filtered_contig_annotations.csv"))
COV_09_19_2023_P3$barcode <- paste0("COV_09_19_2023_P3_", COV_09_19_2023_P3$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_19_2023_P3$contig_id <- paste0("COV_09_19_2023_P3_", COV_09_19_2023_P3$contig_id)

COV_09_19_2023P2P3 <- rbind(COV_09_19_2023_P2, COV_09_19_2023_P3)

#4951574848
indices <- unlist(sapply(s4951574848$barcode, function(x){
  grep(x,COV_09_19_2023P2P3$barcode)
}))
csv_4951574848 <- COV_09_19_2023P2P3[indices,]
write.csv(csv_4951574848, file = here::here("03_Immcantation","s4951574848_filtered_contig_annotations.csv"))

#5357484948
indices <- unlist(sapply(s5357484948$barcode, function(x){
  grep(x,COV_09_19_2023P2P3$barcode)
}))
csv_5357484948 <- COV_09_19_2023P2P3[indices,]
write.csv(csv_5357484948, file = here::here("03_Immcantation","s5357484948_filtered_contig_annotations.csv"))

rm(COV_09_19_2023_P2)
rm(COV_09_19_2023_P3)
rm(COV_09_19_2023P2P3)
#####

#####
# "COV_09_20_2023_P1"
COV_09_20_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_09_20_2023_P1$barcode <- paste0("COV_09_20_2023_P1_", COV_09_20_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_20_2023_P1$contig_id <- paste0("COV_09_20_2023_P1_", COV_09_20_2023_P1$contig_id)

#4954554848
indices <- unlist(sapply(s4954554848$barcode, function(x){
  grep(x,COV_09_20_2023_P1$barcode)
}))
csv_4954554848 <- COV_09_20_2023_P1[indices,]
write.csv(csv_4954554848, file = here::here("03_Immcantation","s4954554848_filtered_contig_annotations.csv"))

#4957484948
indices <- unlist(sapply(s4957484948$barcode, function(x){
  grep(x,COV_09_20_2023_P1$barcode)
}))
csv_4957484948 <- COV_09_20_2023_P1[indices,]
write.csv(csv_4957484948, file = here::here("03_Immcantation","s4957484948_filtered_contig_annotations.csv"))

rm(COV_09_20_2023_P1)
#####

#####
# "COV_09_20_2023_P2/P3"
COV_09_20_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_20_2023_P2$barcode <- paste0("COV_09_20_2023_P2_", COV_09_20_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_20_2023_P2$contig_id <- paste0("COV_09_20_2023_P2_", COV_09_20_2023_P2$contig_id)

COV_09_20_2023_P3 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-20-2023-VDJ-Pool3","outs","filtered_contig_annotations.csv"))
COV_09_20_2023_P3$barcode <- paste0("COV_09_20_2023_P3_", COV_09_20_2023_P3$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_20_2023_P3$contig_id <- paste0("COV_09_20_2023_P3_", COV_09_20_2023_P3$contig_id)

COV_09_20_2023P2P3 <- rbind(COV_09_20_2023_P2, COV_09_20_2023_P3)

#4848544848
indices <- unlist(sapply(s4848544848$barcode, function(x){
  grep(x,COV_09_20_2023P2P3$barcode)
}))
csv_4848544848 <- COV_09_20_2023P2P3[indices,]
write.csv(csv_4848544848, file = here::here("03_Immcantation","s4848544848_filtered_contig_annotations.csv"))

#5053564848 
indices <- unlist(sapply(s5053564848 $barcode, function(x){
  grep(x,COV_09_20_2023P2P3$barcode)
}))
csv_5053564848  <- COV_09_20_2023P2P3[indices,]
write.csv(csv_5053564848, file = here::here("03_Immcantation","s5053564848_filtered_contig_annotations.csv"))

rm(COV_09_20_2023_P2)
rm(COV_09_20_2023_P3)
rm(COV_09_20_2023P2P3)
#####

#####
# "COV_09_21_2023_P1"
COV_09_21_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-21-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_09_21_2023_P1$barcode <- paste0("COV_09_21_2023_P1_", COV_09_21_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_21_2023_P1$contig_id <- paste0("COV_09_21_2023_P1_", COV_09_21_2023_P1$contig_id)

#4955534848
indices <- unlist(sapply(s4955534848$barcode, function(x){
  grep(x,COV_09_21_2023_P1$barcode)
}))
csv_4955534848 <- COV_09_21_2023_P1[indices,]
write.csv(csv_4955534848, file = here::here("03_Immcantation","s4955534848_filtered_contig_annotations.csv"))

#s4957564848
indices <- unlist(sapply(s4957564848$barcode, function(x){
  grep(x,COV_09_21_2023_P1$barcode)
}))
csv_s4957564848 <- COV_09_21_2023_P1[indices,]
write.csv(csv_s4957564848, file = here::here("03_Immcantation","s4957564848_filtered_contig_annotations.csv"))

rm(COV_09_21_2023_P1) #so I don't get confused
#####

#####
# "COV_09_21_2023_P2"
COV_09_21_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-21-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_21_2023_P2$barcode <- paste0("COV_09_21_2023_P2_", COV_09_21_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_21_2023_P2$contig_id <- paste0("COV_09_21_2023_P2_", COV_09_21_2023_P2$contig_id)

#4953494948
indices <- unlist(sapply(s4953494948$barcode, function(x){
  grep(x,COV_09_21_2023_P2$barcode)
}))
csv_4953494948 <- COV_09_21_2023_P2[indices,]
write.csv(csv_4953494948, file = here::here("03_Immcantation","s4953494948_filtered_contig_annotations.csv"))

#5553564848
indices <- unlist(sapply(s5553564848$barcode, function(x){
  grep(x,COV_09_21_2023_P2$barcode)
}))
csv_s5553564848 <- COV_09_21_2023_P2[indices,]
write.csv(csv_s5553564848, file = here::here("03_Immcantation","s5553564848_filtered_contig_annotations.csv"))

rm(COV_09_21_2023_P2)
#####

#####
# "COV_09_26_2023_P1"
COV_09_26_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-26-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_09_26_2023_P1$barcode <- paste0("COV_09_26_2023_P1_", COV_09_26_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_26_2023_P1$contig_id <- paste0("COV_09_26_2023_P1_", COV_09_26_2023_P1$contig_id)

#5150564848
indices <- unlist(sapply(s5150564848$barcode, function(x){
  grep(x,COV_09_26_2023_P1$barcode)
}))
csv_5150564848 <- COV_09_26_2023_P1[indices,]
write.csv(csv_5150564848, file = here::here("03_Immcantation","s5150564848_filtered_contig_annotations.csv"))

#5557544848
indices <- unlist(sapply(s5557544848$barcode, function(x){
  grep(x,COV_09_26_2023_P1$barcode)
}))
csv_s5557544848 <- COV_09_26_2023_P1[indices,]
write.csv(csv_s5557544848, file = here::here("03_Immcantation","s5557544848_filtered_contig_annotations.csv"))

rm(COV_09_26_2023_P1)
#####

#####
# "COV_09_26_2023_P2"
COV_09_26_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-26-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_26_2023_P2$barcode <- paste0("COV_09_26_2023_P2_", COV_09_26_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_26_2023_P2$contig_id <- paste0("COV_09_26_2023_P2_", COV_09_26_2023_P2$contig_id)

#5050484948
indices <- unlist(sapply(s5050484948$barcode, function(x){
  grep(x,COV_09_26_2023_P2$barcode)
}))
csv_5050484948 <- COV_09_26_2023_P2[indices,]
write.csv(csv_5050484948, file = here::here("03_Immcantation","s5050484948_filtered_contig_annotations.csv"))

#5048544848
indices <- unlist(sapply(s5048544848$barcode, function(x){
  grep(x,COV_09_26_2023_P2$barcode)
}))
csv_s5048544848 <- COV_09_26_2023_P2[indices,]
write.csv(csv_s5048544848, file = here::here("03_Immcantation","s5048544848_filtered_contig_annotations.csv"))

rm(COV_09_26_2023_P2)
#####

#####
# "COV_09_27_2023_P1"
COV_09_27_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-27-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_09_27_2023_P1$barcode <- paste0("COV_09_27_2023_P1_", COV_09_27_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_27_2023_P1$contig_id <- paste0("COV_09_27_2023_P1_", COV_09_27_2023_P1$contig_id)

#5750564848
indices <- unlist(sapply(s5750564848$barcode, function(x){
  grep(x,COV_09_27_2023_P1$barcode)
}))
csv_5750564848 <- COV_09_27_2023_P1[indices,]
write.csv(csv_5750564848, file = here::here("03_Immcantation","s5750564848_filtered_contig_annotations.csv"))

#5755544848
indices <- unlist(sapply(s5755544848$barcode, function(x){
  grep(x,COV_09_27_2023_P1$barcode)
}))
csv_s5755544848 <- COV_09_27_2023_P1[indices,]
write.csv(csv_s5755544848, file = here::here("03_Immcantation","s5755544848_filtered_contig_annotations.csv"))

rm(COV_09_27_2023_P1)
#####

#####
# "COV_09_27_2023_P2"
COV_09_27_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-27-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_27_2023_P2$barcode <- paste0("COV_09_27_2023_P2_", COV_09_27_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_27_2023_P2$contig_id <- paste0("COV_09_27_2023_P2_", COV_09_27_2023_P2$contig_id)

#4950544848
indices <- unlist(sapply(s4950544848$barcode, function(x){
  grep(x,COV_09_27_2023_P2$barcode)
}))
csv_4950544848 <- COV_09_27_2023_P2[indices,]
write.csv(csv_4950544848, file = here::here("03_Immcantation","s4950544848_filtered_contig_annotations.csv"))

#5249544848
indices <- unlist(sapply(s5249544848$barcode, function(x){
  grep(x,COV_09_27_2023_P2$barcode)
}))
csv_s5249544848 <- COV_09_27_2023_P2[indices,]
write.csv(csv_s5249544848, file = here::here("03_Immcantation","s5249544848_filtered_contig_annotations.csv"))

rm(COV_09_27_2023_P2)
#####

#####
# "COV_09_28_2023_P1/P3"
COV_09_28_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_09_28_2023_P1$barcode <- paste0("COV_09_28_2023_P1_", COV_09_28_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_28_2023_P1$contig_id <- paste0("COV_09_28_2023_P1_", COV_09_28_2023_P1$contig_id)

COV_09_28_2023_P3 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool3","outs","filtered_contig_annotations.csv"))
COV_09_28_2023_P3$barcode <- paste0("COV_09_28_2023_P3_", COV_09_28_2023_P3$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_28_2023_P3$contig_id <- paste0("COV_09_28_2023_P3_", COV_09_28_2023_P3$contig_id)

COV_09_28_2023P1P3 <- rbind(COV_09_28_2023_P1, COV_09_28_2023_P3)

#5054574848
indices <- unlist(sapply(s5054574848$barcode, function(x){
  grep(x,COV_09_28_2023P1P3$barcode)
}))
csv_5054574848 <- COV_09_28_2023P1P3[indices,]
write.csv(csv_5054574848, file = here::here("03_Immcantation","s5054574848_filtered_contig_annotations.csv"))

#5456544848
indices <- unlist(sapply(s5456544848$barcode, function(x){
  grep(x,COV_09_28_2023P1P3$barcode)
}))
csv_5456544848 <- COV_09_28_2023P1P3[indices,]
write.csv(csv_5456544848, file = here::here("03_Immcantation","s5456544848_filtered_contig_annotations.csv"))

rm(COV_09_28_2023_P1)
rm(COV_09_28_2023_P3)
rm(COV_09_28_2023P2P3)
#####

#####
# "COV_09_28_2023_P2"
COV_09_28_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-09-28-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_09_28_2023_P2$barcode <- paste0("COV_09_28_2023_P2_", COV_09_28_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_09_28_2023_P2$contig_id <- paste0("COV_09_28_2023_P2_", COV_09_28_2023_P2$contig_id)

#5356534848
indices <- unlist(sapply(s5356534848$barcode, function(x){
  grep(x,COV_09_28_2023_P2$barcode)
}))
csv_5356534848 <- COV_09_28_2023_P2[indices,]
write.csv(csv_5356534848, file = here::here("03_Immcantation","s5356534848_filtered_contig_annotations.csv"))

#5556484948
indices <- unlist(sapply(s5556484948$barcode, function(x){
  grep(x,COV_09_28_2023_P2$barcode)
}))
csv_s5556484948 <- COV_09_28_2023_P2[indices,]
write.csv(csv_s5556484948, file = here::here("03_Immcantation","s5556484948_filtered_contig_annotations.csv"))

rm(COV_09_28_2023_P2)
#####

#####
# "COV_10_05_2023_P1"
COV_10_05_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-10-05-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_10_05_2023_P1$barcode <- paste0("COV_10_05_2023_P1_", COV_10_05_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_10_05_2023_P1$contig_id <- paste0("COV_10_05_2023_P1_", COV_10_05_2023_P1$contig_id)

#5553554848
indices <- unlist(sapply(s5553554848$barcode, function(x){
  grep(x,COV_10_05_2023_P1$barcode)
}))
csv_5553554848 <- COV_10_05_2023_P1[indices,]
write.csv(csv_5553554848, file = here::here("03_Immcantation","s5553554848_filtered_contig_annotations.csv"))

rm(COV_10_05_2023_P1)
#####

#####
# "COV_10_06_2023_P1"
COV_10_06_2023_P1 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-10-06-2023-VDJ-Pool1","outs","filtered_contig_annotations.csv"))
COV_10_06_2023_P1$barcode <- paste0("COV_10_06_2023_P1_", COV_10_06_2023_P1$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_10_06_2023_P1$contig_id <- paste0("COV_10_06_2023_P1_", COV_10_06_2023_P1$contig_id)

#5657574848
indices <- unlist(sapply(s5657574848$barcode, function(x){
  grep(x,COV_10_06_2023_P1$barcode)
}))
csv_5657574848 <- COV_10_06_2023_P1[indices,]
write.csv(csv_5657574848, file = here::here("03_Immcantation","s5657574848_filtered_contig_annotations.csv"))

rm(COV_10_06_2023_P1)
#####

#####
# "COV_10_06_2023_P2"
COV_10_06_2023_P2 <- read.csv(file = here::here("02_CellRanger","VRC-COVAIL-10-06-2023-VDJ-Pool2","outs","filtered_contig_annotations.csv"))
COV_10_06_2023_P2$barcode <- paste0("COV_10_06_2023_P2_", COV_10_06_2023_P2$barcode) #go ahead and paste sequence name so we don't need the 01 script
COV_10_06_2023_P2$contig_id <- paste0("COV_10_06_2023_P2_", COV_10_06_2023_P2$contig_id)


#4856554848
indices <- unlist(sapply(s4856554848$barcode, function(x){
  grep(x,COV_10_06_2023_P2$barcode)
}))
csv_4856554848 <- COV_10_06_2023_P2[indices,]
write.csv(csv_4856554848, file = here::here("03_Immcantation","s4856554848_filtered_contig_annotations.csv"))

#5653544848
indices <- unlist(sapply(s5653544848$barcode, function(x){
  grep(x,COV_10_06_2023_P2$barcode)
}))
csv_s5653544848 <- COV_10_06_2023_P2[indices,]
write.csv(csv_s5653544848, file = here::here("03_Immcantation","s5653544848_filtered_contig_annotations.csv"))

rm(COV_10_06_2023_P2)
#####

# paths <- c("../03_Immcantation/s5457484948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5048574848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4951574848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4954554848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5357484948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4954554848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4957484948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4848544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5053564848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4955534848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4957564848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4953494948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5553564848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5150564848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5557544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5050484948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5048544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5750564848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5755544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4950544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5249544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5054574848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5456544848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5356534848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5556484948_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5553554848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5657574848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s4856554848_filtered_contig_annotations.csv",
#            "../03_Immcantation/s5653544848_filtered_contig_annotations.csv"
# )

sample_names <- c("s5351564848",
                "s5457484948",
                "s5048574848",
                "s4951574848",
                "s4954554848",
                "s5357484948",
                "s4954554848",
                "s4957484948",
                "s4848544848",
                "s5053564848",
                "s4955534848",
                "s4957564848",
                "s4953494948",
                "s5553564848",
                "s5150564848",
                "s5557544848",
                "s5050484948",
                "s5048544848",
                "s5750564848",
                "s5755544848",
                "s4950544848",
                "s5249544848",
                "s5054574848",
                "s5456544848",
                "s5356534848",
                "s5556484948",
                "s5553554848",
                "s5657574848",
                "s4856554848",
                "s5653544848"
)

# for(i in 1:length(paths)){
#   x <- read.csv(file = paths[i])
#   x$barcode <- paste0(sample_names[i],"_",x$barcode)
#   x$contig_id <- paste0(sample_names[i],"_",x$contig_id)
#   write.csv(x, file = paste0(sample_names[i],"_filtered_contig_annotations.csv"))
# }

fasta <- paste0(sample_names,"_filtered_contig.fasta")
annotations <- paste0(sample_names,"_filtered_contig_annotations.csv")
samp_tab <- data.frame(sample_names,fasta,annotations)

write.table(samp_tab, file = "vdj_files.txt", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote=FALSE)














# concatenate A and B so that there's only 2 files per sample (M and P)
# remove MemA and MemB files
# M415A <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S8_316_415_MemA_VDJLib/outs/M415A_filtered_contig_annotations.csv")
# M415A$barcode <- paste0("M415A","_",M415A$barcode)
# M415A$contig_id <- paste0("M415A","_",M415A$contig_id)
# 
# M415B <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S9_316_415_MemB_VDJLib/outs/M415B_filtered_contig_annotations.csv")
# M415B$barcode <- paste0("M415B","_",M415B$barcode)
# M415B$contig_id <- paste0("M415B","_",M415B$contig_id)
# 
# M415 <- rbind(M415A,M415B)
# 
# ## add in plasmablasts
# P415 <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S10_316_415_PB_VDJLib/outs/P415_filtered_contig_annotations.csv")
# P415$barcode <- paste0("P415","_",P415$barcode)
# P415$contig_id <- paste0("P415","_",P415$contig_id)
# MP415 <- rbind(M415,P415)
# write.csv(MP415, file = "MP415_filtered_contig_annotations.csv")
# 
# 
# M402A <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S5_316_402_MemA_VDJLib/outs/M402A_filtered_contig_annotations.csv")
# M402A$barcode <- paste0("M402A","_",M402A$barcode)
# M402A$contig_id <- paste0("M402A","_",M402A$contig_id)
# 
# M402B <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S6_316_402_MemB_VDJLib/outs/M402B_filtered_contig_annotations.csv")
# M402B$barcode <- paste0("M402B","_",M402B$barcode)
# M402B$contig_id <- paste0("M402B","_",M402B$contig_id)
# 
# M402 <- rbind(M402A,M402B)
# 
# P402 <- read.csv(file = "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S7_316_402_PB_VDJLib/outs/P402_filtered_contig_annotations.csv")
# P402$barcode <- paste0("P402","_",P402$barcode)
# P402$contig_id <- paste0("P402","_",P402$contig_id)
# MP402 <- rbind(M402,P402)
# 
# write.csv(MP402, file = "MP402_filtered_contig_annotations.csv")
# 
# 
#     m_paths <- c("../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S3_316_201_Mem_VDJLib/outs/M201_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S1_316_416_Mem_VDJLib/outs/M416_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S11_316_418_Mem_VDJLib/outs/M418_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S16_316_410_Mem_5_VDJLib/outs/filtered_contig_annotations.csv",
# "../../ver_020/S1_A1_408_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S3_A3_412_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S5_A5_403_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S10_A10_405_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S11_A11_309_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S13_B1_302_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv",
#            "../../ver_020/S15_B3_318_Mem_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv")
# 
# p_paths <- c("../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S4_316_201_PB_VDJLib/outs/P201_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S2_316_416_PB_VDJLib/outs/P416_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S12_316_418_PB_VDJLib/outs/P418_filtered_contig_annotations.csv",
# "../../10XNxtSq_CITESeq_Exp19_210825_BH3CGHDMXY/ver_010_LATEST_Transfer/S14_316_410_PB_2_VDJLib/outs/P410_filtered_contig_annotations.csv",
# "../../ver_020/S2_A2_408_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
# "../../ver_020/S4_A4_412_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
# "../../ver_020/S6_A6_403_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
#  "../../ver_020/S9_A9_405_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
#  "../../ver_020/S12_A12_309_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
#  "../../ver_020/S14_B2_302_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv",
#  "../../ver_020/S16_B4_318_PB_Plasmablasts_VDJLib/outs/filtered_contig_annotations.csv")
# 
# subjects<- c("201","416","418","410","408","412","403","405","309","302","318")
# 
# samp_tab <- data.frame(m_paths,p_paths,subjects)
# 
# for(i in 1:nrow(samp_tab)){
#   mem_ann <- read.csv(file = samp_tab[i,"m_paths"])
#   mem_ann$barcode <- paste0("M",samp_tab[i,"subjects"],"_",mem_ann$barcode)
#   mem_ann$contig_id <- paste0("M",samp_tab[i,"subjects"],"_",mem_ann$contig_id)
# 
#   p_ann <- read.csv(file = samp_tab[i,"p_paths"])
#   p_ann$barcode <- paste0("P",samp_tab[i,"subjects"],"_",p_ann$barcode)
#   p_ann$contig_id <- paste0("P",samp_tab[i,"subjects"],"_",p_ann$contig_id)
#   
#   x <- rbind(mem_ann,p_ann)
#   write.csv(x, file = paste0("MP",samp_tab[i,"subjects"],"_filtered_contig_annotations.csv"))
# }
# 
# # now process memA and memB pools and concatenate them to their corresponding files
# memA <- read.csv(file = "../../ver_020/S7_A7_Mem_A_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv")
# memB <- read.csv(file = "../../ver_020/S8_A8_Mem_B_Mem_B_cells_VDJLib/outs/filtered_contig_annotations.csv")
# 
# memA$barcode <- paste0("MemA","_",memA$barcode)
# memA$contig_id <- paste0("MemA","_",memA$contig_id)
# 
# memB$barcode <- paste0("MemB","_",memB$barcode)
# memB$contig_id <- paste0("MemB","_",memB$contig_id)
# 
# Mem <- rbind(memA,memB)
# 
# ## read in barcode files for each subejct 
# S410 <- read.csv(file = "MemA_410_barcodes.csv", row.names = 1)
# S416 <- read.csv(file = "MemA_416_barcodes.csv", row.names = 1)
# S418 <- read.csv(file = "MemA_418_barcodes.csv", row.names = 1)
# S201 <- read.csv(file = "MemB_201_barcodes.csv", row.names = 1)
# S402 <- read.csv(file = "MemB_402_barcodes.csv", row.names = 1)
# S415 <- read.csv(file = "MemB_415_barcodes.csv", row.names = 1)
# 
# #410
# indices <- unlist(sapply(S410$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_410 <- Mem[indices,]
# Orig_410 <- read.csv(file = "MP410_filtered_contig_annotations.csv", row.names=1)
# New_410 <- rbind(Mem_410,Orig_410)
# write.csv(New_410, file = "MP410_filtered_contig_annotations.csv")
# 
# #416
# indices <- unlist(sapply(S416$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_416 <- Mem[indices,]
# Orig_416 <- read.csv(file = "MP416_filtered_contig_annotations.csv", row.names=1)
# New_416 <- rbind(Mem_416,Orig_416)
# write.csv(New_416, file = "MP416_filtered_contig_annotations.csv")
# 
# #418
# indices <- unlist(sapply(S418$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_418 <- Mem[indices,]
# Orig_418 <- read.csv(file = "MP418_filtered_contig_annotations.csv", row.names=1)
# New_418 <- rbind(Mem_418,Orig_418)
# write.csv(New_418, file = "MP418_filtered_contig_annotations.csv")
# 
# 
# #201
# indices <- unlist(sapply(S201$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_201 <- Mem[indices,]
# Orig_201 <- read.csv(file = "MP201_filtered_contig_annotations.csv", row.names=1)
# New_201 <- rbind(Mem_201,Orig_201)
# write.csv(New_201, file = "MP201_filtered_contig_annotations.csv")
# 
# #402
# indices <- unlist(sapply(S402$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_402 <- Mem[indices,]
# Orig_402 <- read.csv(file = "MP402_filtered_contig_annotations.csv", row.names=1)
# New_402 <- rbind(Mem_402,Orig_402)
# write.csv(New_402, file = "MP402_filtered_contig_annotations.csv")
# 
# #415
# indices <- unlist(sapply(S415$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_415 <- Mem[indices,]
# Orig_415 <- read.csv(file = "MP415_filtered_contig_annotations.csv", row.names=1)
# New_415 <- rbind(Mem_415,Orig_415)
# write.csv(New_415, file = "MP415_filtered_contig_annotations.csv")
# 
# ##run5 M309B 
# ann_309 <- read.csv(file = here::here("analysis","VDJ_run3_run4_by_subject","MP309_filtered_contig_annotations.csv"), row.names = 1)
# # ann_309 <- ann_309[,2:32] #remove X column
# extra_ann_309 <- read.csv(file = here::here("CellRanger_run5","S-316-309-MBC-VDJ","outs","filtered_contig_annotations.csv"))
# extra_ann_309$barcode <- paste0("M309B_",extra_ann_309$barcode)
# extra_ann_309$contig_id <- paste0("M309B_",extra_ann_309$contig_id)
# 
# ann_309 <- rbind(ann_309,extra_ann_309)
# write.csv(ann_309,file = here::here("analysis","VDJ_run3_run4_by_subject","MP309_filtered_contig_annotations.csv"))
# 
# ##split run5 M302318 pool annotation files and add them to other annotation files
# Mem <- read.csv(file = here::here("CellRanger_run5","S-316-302-318-Mem-VDJ","outs","filtered_contig_annotations.csv"))
# 
# Mem$barcode <- paste0("M302318","_",Mem$barcode)
# Mem$contig_id <- paste0("M302318","_",Mem$contig_id)
# 
# ## read in barcode files for each subejct 
# S302 <- read.csv(file = here::here("analysis","VDJ_run3_run4_by_subject","M302318_302_barcodes.csv"))
# S318 <- read.csv(file = here::here("analysis","VDJ_run3_run4_by_subject","M302318_318_barcodes.csv"))
# 
# #302
# indices <- unlist(sapply(S302$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_302 <- Mem[indices,]
# Orig_302 <- read.csv(file = here::here("analysis","VDJ_run3_run4_by_subject","MP302_filtered_contig_annotations.csv"), row.names=1)
# New_302 <- rbind(Mem_302,Orig_302)
# write.csv(New_302, file = here::here("analysis","VDJ_run3_run4_by_subject","MP302_filtered_contig_annotations.csv"))
# 
# #318
# indices <- unlist(sapply(S318$barcode, function(x){
#     grep(x,Mem$barcode)
# }))
# Mem_318 <- Mem[indices,]
# Orig_318 <- read.csv(file = here::here("analysis","VDJ_run3_run4_by_subject","MP318_filtered_contig_annotations.csv"), row.names=1)
# New_318 <- rbind(Mem_318,Orig_318)
# write.csv(New_318, file = here::here("analysis","VDJ_run3_run4_by_subject","MP318_filtered_contig_annotations.csv"))
# 
# sample_names <- c("MP302","MP309","MP318")
# fasta <- c("MP302_filtered_contig.fasta","MP309_filtered_contig.fasta","MP318_filtered_contig.fasta")
# annotations <- c("MP302_filtered_contig_annotations.csv","MP309_filtered_contig_annotations.csv","MP318_filtered_contig_annotations.csv")
# samp_tab <- data.frame(sample_names,fasta,annotations)
# 
# write.table(samp_tab, file = here::here("analysis","VDJ_run3_run4_by_subject","vdj_files_run5.txt"), sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote=FALSE)