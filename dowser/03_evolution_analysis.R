# Have to load Apptainer container on Skyline then run R
library(dowser)
library(dplyr)
library(ggtree)
library(sessioninfo)

#load trees with meta data
trees <- readRDS(file = here::here("dowser","HL_trees_Timepoint_raxml.rds"))
results <- readRDS(file = here::here("dowser","createGermlines.rds"))

### measuring evolution ###
evo.trees <- correlationTest(trees,time = "Timepoint.num",permutations = 10000,perm_type = "uniform")
evo.trees <- evo.trees[order(evo.trees$p),]

### add meta data to trees object results$clone_subgroup_id = trees$clone_id
meta <- results %>% select(clone_subgroup_id, Booster, adj.ProtoOmi,Infection) %>% distinct(clone_subgroup_id, .keep_all = TRUE)
meta$clone_id <- meta$clone_subgroup_id
meta <- meta %>% filter(clone_id %in% evo.trees$clone_id)
evo.trees <- left_join(evo.trees,meta, by = "clone_id")
saveRDS(evo.trees,file = here::here("dowser","HL_trees_Timepoint_raxml_Evo.rds"))
# evo.trees <- readRDS(file = here::here("dowser","HL_trees_Timepoint_raxml_Evo.rds"))

## make csv of lineages
dat <- evo.trees %>% select(clone_id,locus,seqs,slope,correlation,p,Booster,adj.ProtoOmi,Booster,Infection)
dat <- dat %>% mutate(sig = case_when(p < 0.05 ~ TRUE,
                              p > 0.05 ~ FALSE))
write.csv(dat,file = here::here("dowser","Evo_dat_Timepoint_uniform.csv"),row.names = FALSE)
# dat <- read.csv(file = here::here("dowser","Evo_dat_Timepoint_uniform.csv"))

sig_counts <- dat %>% filter(sig == TRUE) %>% count(Booster,Infection)
# # A tibble: 6 × 3
#   Booster               Infection     n
#   <chr>                 <chr>     <int>
# 1 Omicron               N             3
# 2 Omicron               Y             4
# 3 Omicron And Prototype N             5
# 4 Omicron And Prototype Y             2
# 5 Prototype             N             5
# 6 Prototype             Y             4

# > table(dat$Booster)

#               Omicron Omicron And Prototype             Prototype 
#                   266                   208                   235 

# > 7/266
# [1]0.02631579
# > 7/208
# [1] 0.03365385
# > 9/235
# [1] 0.03829787

min_seq <- dat %>% filter(sig == TRUE) %>% summarise(min_seq = min(seqs))
print(min_seq)
# 4

pdf(file = here::here("dowser","Evo_scatterplot_uniform.pdf"))
ggplot(dat, aes(x = Booster, y = slope, colour = sig, shape = adj.ProtoOmi)) + 
geom_point(position = position_jitter(), size = 1) +
scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "gray"), name = "Significant") +
scale_shape_discrete(name = "Specificity") +
theme_bw() +
theme(axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(dat, aes(x = Booster, y = correlation, colour = sig, shape = adj.ProtoOmi)) + 
geom_point(position = position_jitter(), size = 1) +
scale_colour_manual(values = c("TRUE" = "black", "FALSE" = "gray"), name = "Significant") +
scale_shape_discrete(name = "Specificity") +
theme_bw() +
theme(axis.text = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

##### 2 ##### Timepoint.num
plots <- plotTrees(trees, tips="Timepoint.num",tip_palette = "RdYlBu", tipsize = 3)
plots <- lapply(plots, function(x){
    x + geom_tiplab(aes(label = Timepoint),offset = 0.03)
})

treesToPDF(plots, here::here("dowser","tree_plots_raxml_timepoint.pdf"), nrow = 1, ncol = 2, height = 20)


## plots for manuscript ##
# s5048544848_374_11_1 evolving
# s5553564848_1439_26_1 not evolving 
idx <- c(grep("s5048544848_374_11_1", trees$clone_id),grep("s5553564848_1439_26_1", trees$clone_id))
treesToPDF(plots[idx], here::here("dowser","tree_plots_raxml_timepoint_manuscript.pdf"), nrow = 1, ncol = 1, height = 10)


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()