library(Seurat)
library(tidyverse)
library(sessioninfo)
library(tidyseurat)
library(here)

mergedPosObject.dsb <-readRDS(file = here::here("04_Analysis","data_objects","02_dsb_normalization","MergedSeuratObject_positive_dsbnormalized.rds"))
####################################### Look at DSB Normalized Protein Expression################################

#set run num and load in csv
runNum <- 1
run_info <- read.csv(file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv")) #still having issues finding out how to make this sheet from Abby's pipeline

mergedPosObject.dsb %>% count(Timepoint)

pdf(file = here::here("04_Analysis","plots","02_dsb_normalization","CD14_v_CD19_all.pdf"))
mergedPosObject.dsb %>%
  join_features(features = c("P-CD14","P-CD19")) %>% 
  select(one_of(c(".cell", ".feature",".abundance_Prot","Timepoint"))) %>% 
  pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
  ggplot(aes(x= `P-CD14`, y=`P-CD19`)) + 
  geom_point(aes(color = Timepoint), size=0.2) +
  scale_x_continuous(limits = c(-2, 25)) +
  scale_y_continuous(limits = c(-2, 25)) +
  facet_wrap(.~Timepoint) +
  theme(aspect.ratio=1) +
  theme_bw()
dev.off()


#Label Bcells in Meta data
mergedPosObject.dsb$Bcell <- mergedPosObject.dsb %>% 
  tidyseurat::join_features(features = c("P-CD19","P-CD3","P-CD14","P-CD56")) %>% 
  pivot_wider(names_from = .feature, values_from = .abundance_Prot) %>% 
  mutate(Bcell = if_else(`P-CD19` > 2.5 & `P-CD3` < 7.5 & `P-CD14` < 5 & `P-CD56` < 10, TRUE, FALSE)) %>%
  pull("Bcell")

run_info$bcells[runNum] <- table(mergedPosObject.dsb$Bcell)[2]
message("Number of B cells total (TRUE indicates B-Cell status): \n",paste0(capture.output(table(mergedPosObject.dsb$Bcell)), collapse = "\n"))
message("Number of B cells per pool: \n",paste0(capture.output(table(mergedPosObject.dsb$Bcell, mergedPosObject.dsb$orig.ident)), collapse = "\n"))

write.csv(run_info,file = here::here("04_Analysis","data_objects","01_build_seurat","run_info_stats.csv"),row.names=FALSE)
saveRDS(mergedPosObject.dsb, file = here::here("04_Analysis","data_objects","02_dsb_normalization","MergedSeuratObject_positive_dsbnormalized_BCellAnnotated.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()