library(phylotools)
library(dplyr)
library(readxl)
library(here)
library(Seurat)
library(tidyr)

#read in object
seuObj <- readRDS(file = here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "covObj_clustered_demultiplexed.rds"))
df <- seuObj@meta.data %>% select(., CELL, LC_locus, sequence, LC_sequence)

#read in excel files with antibodies of choice
omicronAbs <- read_xlsx(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "AntibodiesChosenSoFar_better_donor_coverage.xlsx"), sheet = "pickedomicron")
protoSpecific <- read_xlsx(here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "AntibodiesChosenSoFar_better_donor_coverage.xlsx"), sheet = "pickedOmiActivated")

#convert to long format
fasta <- df %>% filter(CELL %in% c(omicronAbs$CELL, protoSpecific$CELL)) %>%
      pivot_longer(!c(CELL, LC_locus), names_to = "sequence") %>%
      mutate(CELL = case_when(sequence == "sequence" ~ paste(CELL, "HC", sep = "_"),
             TRUE ~ paste(CELL, LC_locus, sep="_"))) %>%
      rename(`seq.name` = CELL, `seq.text` = value)

#write to a fasta file
#we need to split in two because imgt only wants 50 at a time :(
dat2fasta(fasta[fasta$sequence == "sequence", colnames(fasta) %in% c("seq.name", "seq.text")], here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "antibody_fastas_HC.fasta"))
dat2fasta(fasta[fasta$sequence == "LC_sequence", colnames(fasta) %in% c("seq.name", "seq.text")], here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "antibody_fastas_LC.fasta"))
