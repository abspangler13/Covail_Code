library(dowser)
library(dplyr)
library(sessioninfo)

trees <- readRDS(file = here::here("dowser","HL_trees_Specificity_raxml.rds"))
results <- readRDS(file = here::here("dowser","createGermlines.rds"))

### add meta data to trees object
results$clone_id_orig <- results$clone_id
results$clone_id <- results$clone_subgroup_id

## add in other meta data from results object
results <- results %>% 
    select(Booster, ProtoOmi, Infection, clone_id) %>% 
    distinct(clone_id, .keep_all = TRUE)

trees <- left_join(trees, results, by = "clone_id")
saveRDS(trees, file = here::here("dowser","HL_trees_Specificity_raxml_meta.rds"))
# trees <- readRDS(file = here::here("dowser","HL_trees_Cluster_raxml_meta.rds"))


trait = "ProtoOmi"
igphyml_location = "/usr/local/share/igphyml/src/igphyml"

###run find switches on everything together 
switches <- findSwitches(trees,permutations=100, trait=trait, 
  igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)

saveRDS(switches, file = here::here("dowser","switches_specificity_all.rds"))

# ### find switches for cluster
# trees <- readRDS(file = here::here("dowser","HL_trees_Cluster_raxml.rds"))
# trees <- left_join(trees, results, by = "clone_id")
# saveRDS(trees, file = here::here("dowser","HL_trees_Cluster_raxml_meta.rds"))

# trait = "ClusterLabel.AS"
# igphyml_location = "/usr/local/share/igphyml/src/igphyml"

# ###run find switches on everything together 
# switches <- findSwitches(trees,permutations=100, trait=trait, 
#   igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)

# saveRDS(switches, file = here::here("dowser","switches_cluster_all.rds"))

# ### divide trees up by Booster
# trees.ls <- split(trees,f = trees$Booster)

# switches.ls <- lapply(trees.ls,function(x){
#   findSwitches(x,permutations=100, trait=trait, 
#   igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
# })

# saveRDS(switches.ls, file = here::here("dowser","switches_cluster_booster.rds"))

# ### divide trees up by broad specificity AND exposure
# trees$specificity_exposure <- paste0(trees$Spec.Broad,"_",trees$exposure)
# trees.ls <- split(trees,f = c(trees$specificity_exposure))

# switches.ls <- lapply(trees.ls,function(x){
#   findSwitches(x,permutations=100, trait=trait, 
#   igphyml=igphyml_location, fixtrees=TRUE,force_resolve = TRUE)
# })

# saveRDS(switches.ls, file = here::here("dowser","switches_cluster_specificity_exposure.rds"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()