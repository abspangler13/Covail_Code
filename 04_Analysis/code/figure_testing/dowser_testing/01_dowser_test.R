#demultiplexing infected cohort and labelling phenotypic classifications from Sarah's paper

#load the dependencies
library(ggplot2)
library(dplyr)
library(here)
library(dowser)
library(vroom)
library(stringr)
library(readxl)
library(sessioninfo)

#load the data- we need this as it contains our timepoint info
metadata <- rbind(read_xlsx(path=here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "MonoclonalsToChoose_Uninfected_Infected_PresentAtDay180AndEarlies.xlsx"), sheet="Infected"),
                  read_xlsx(path=here::here("04_Analysis", "data_objects", "06_repertoire_analysis", "MonoclonalsToChoose_Uninfected_Infected_PresentAtDay180AndEarlies.xlsx"), sheet="Infected"))

#read in the vdj data so that we can rebuild the trees
clones <- vroom(paste0(here::here("03_Immcantation"), "/",list.files(path = here::here("03_Immcantation"),pattern = "*heavy_light_clone-pass_germ-pass.tsv")))

#add timepoint to vdj data and remove anything that doesn't have a timepoint
clones$Timepoint <- metadata$Timepoint[match(clones$sequence_id, metadata$sequence_id)]
clones <- clones %>% filter(sequence_id %in% metadata$sequence_id) %>%
            mutate(Timepoint = as.numeric(str_extract(Timepoint, "(?<=Day )[0-9]+")))

#now we can do dowser
clonesFormat = formatClones(clones, traits="Timepoint", columns="Timepoint")
trees = getTrees(clonesFormat, build="pml", nproc=1)

#plot the trees
plots = plotTrees(trees,tips="Timepoint")

pdf(file = here::here("04_Analysis", "plots", "07_dowser", "clonal_trees.pdf"))
for(i in 1:length(plots)){
plots[[i]]
}
dev.off()

#do date randomization test
test = correlationTest(trees, permutations=10000, time="Timepoint")
utest = correlationTest(trees, permutation = 10000, time="Timepoint", perm_type="uniform")
utestDF <- as.data.frame(utest)
write.csv(utestDF, file=here::here("04_Analysis", "data_objects", "07_dowser", "utest_correlations.csv"))

#Reproducibility
session_info()