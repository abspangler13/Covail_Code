#I have to add SHM in a separate script because version control issues prevent me
#from using Seurat and Shazam at the same time >:(
library(tidyverse)
library(readr)
library(plyr)
library(shazam)
library(sessioninfo)
library(vroom)
library(here)

#read in VDJ CSV
All.vdj.M <- read.csv(file = here::here("04_Analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))

#function from shazam. calculates somatic hyper- mutation of vdj. compares sequence of vdj to germline. somatic mutations in vdj. naive bcell has undergone recombination but not somatic hypermutation
All.vdj.M <- observedMutations(All.vdj.M, sequenceColumn = "sequence_alignment",
                                   germlineColumn = "germline_alignment_d_mask",
                                   regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = TRUE)

count <- observedMutations(All.vdj.M, sequenceColumn = "sequence_alignment",
                               germlineColumn = "germline_alignment_d_mask",
                               regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = FALSE)

All.vdj.M$mu_count <- count$mu_count[match(All.vdj.M$sequence_id, count$sequence_id)]

count <- observedMutations(All.vdj.M, sequenceColumn = "LC_sequence_alignment",
                            germlineColumn = "LC_germline_alignment", #I'll need to think about how we can run CreateGermlines on LC data- can we do it on the light_parse files?
                            regionDefinition = IMGT_V_BY_REGIONS, combine = TRUE, frequency = FALSE)

All.vdj.M$LC_mu_count <- count$mu_count[match(All.vdj.M$sequence_id, count$sequence_id)]

rm(count)

#write the file
write.csv(All.vdj.M, file = here::here("04_Analysis","data_objects","03_vdj","All.vdj.sum.M_all.csv"))