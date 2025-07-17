library(Seurat)
library(tidyseurat)
library(tidyverse)

A316 <- readRDS(file = here::here("analysis","data_objects","05_clustering","A316_final.rds")) 

# "v_call" # can pick the first one 
# "LC_v_call" # can pick the first one 
# "c_call" #isotype, chose most prevalent
# "LC_c_call" #shouldn't vary across a clone. should either be kappa or lambda
# "specificity.final" 
# "junction_aa"
# "LC_junction_aa" 

#for each clone_subject_id, count how many of each v_call and then pick the most prevalent v_call
tab <- A316 %>% group_by(clone_subject_id) %>% count(v_call) %>% slice(which.min(n)) %>% select(v_call)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(LC_v_call) %>% slice(which.min(n)) %>% select(LC_v_call)
tab <- left_join(tab,tab.2)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(c_call) %>% slice(which.min(n)) %>% select(c_call)
tab <- left_join(tab,tab.2)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(LC_c_call) %>% slice(which.min(n)) %>% select (LC_c_call)
tab <- left_join(tab,tab.2)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(specificity.final) %>% select(specificity.final)
tab <- left_join(tab,tab.2)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(junction_aa) %>% slice(which.min(n)) %>% select (junction_aa)
tab <- left_join(tab,tab.2)
tab.2 <- A316 %>% group_by(clone_subject_id) %>% count(LC_junction_aa) %>% slice(which.min(n)) %>% select (LC_junction_aa)
tab <- left_join(tab,tab.2)
write.csv(tab, file = here::here("analysis","data_objects","04_HA","clone_subject_id_table.csv"))


# j_aa <- A316 %>% group_by(clone_subject_id) %>% count(junction_aa) %>% arrange(clone_subject_id)
# j_aa <- left_join(j_aa,tab, by = "clone_subject_id")
# write.csv(j_aa, file = here::here("analysis","data_objects","05_clustering","clone_subject_id_j_aa.csv"))

# lc_j_aa <- A316 %>% group_by(clone_subject_id) %>% count(LC_junction_aa) %>% arrange(clone_subject_id)
# lc_j_aa <- left_join(lc_j_aa,tab, by = "clone_subject_id")
# write.csv(lc_j_aa, file = here::here("analysis","data_objects","05_clustering","clone_subject_id_lc_j_aa.csv"))
