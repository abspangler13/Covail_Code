library(tidyverse)
library(dplyr)

old.dat <- read.csv(file = here::here("analysis","data_objects","03_vdj","A316_former.csv"))

old.dat <- old.dat %>% select(CELL,clone_subject_id)
colnames(old.dat) <- c("CELL","old_clone_subject_id")

new.dat <- readRDS(file = here::here("analysis","data_objects","04_HA","A316_specificity.rds"))

new.dat <- as.data.frame(as.matrix(new.dat@meta.data))

lookup <- left_join(new.dat,old.dat, by = "CELL")

write.csv(lookup, file = here::here("analysis","data_objects","04_HA","clone_id_lookup_tab.csv"))

