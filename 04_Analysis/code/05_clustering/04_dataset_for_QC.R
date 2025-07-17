library(Seurat)
library(tidyseurat)
library(tidyverse)
library(ggridges)
library(here)
library(writexl)

#load in the two seurat objects
seurat04res <- readRDS(here::here("04_Analysis", "data_objects", "05_clustering", "res_04", "COVAIL_ReclusteredAzimuth_04res.rds"))
seurat03res <- readRDS(here::here("04_Analysis", "data_objects", "05_clustering", "res_03", "COVAIL_ReclusteredAzimuth_03res.rds"))

#add in the protein data to each
seurat04resDF <- seurat04res %>% tidyseurat::join_features(all = TRUE, assay = "Prot") %>% pivot_wider(names_from=.feature,values_from=.abundance_Prot)
seurat03resDF <- seurat03res %>% tidyseurat::join_features(all = TRUE, assay = "Prot") %>% pivot_wider(names_from=.feature,values_from=.abundance_Prot)

#choose the columns to send to sarah
columns <- c("CELL", "Subject", "Booster", "Timepoint", "Infection", "adj.ProtoOmi", "c_call",
             row.names(seurat03res@assays$Prot@counts))

#pull out idents
cluster3 <- seurat03res@active.ident
cluster4 <- seurat04res@active.ident

#make final object
finalSeurat <- seurat03resDF  %>%
                select(all_of(columns)) %>%
                mutate(res03Cluster = cluster3[match(names(cluster3), .$CELL)],
                       res04Cluster = cluster4[match(names(cluster4), .$CELL)])

#write
write_xlsx(finalSeurat, here::here("04_Analysis", "data_objects", "05_clustering", "COVAILSeuratObject_Metadata_Protein_ClusterLabels.xlsx"))

#####make bar graph
bargraph03 <- finalSeurat %>%
              group_by(Timepoint, adj.ProtoOmi, res03Cluster) %>%
              summarize(n = n()) %>%
              mutate(Proportion = n / sum(n),
                     Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))

ggplot(bargraph03, aes(fill = res03Cluster, y = Proportion, x = Timepoint))+
  geom_bar(position = "stack", stat = "identity")+
  facet_grid(cols = vars(adj.ProtoOmi))+
  ggtitle("Resolution 0.3")+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(here::here("04_Analysis", "plots", "05_clustering", "Resolution03_ClustersOverTime.png"), width = 5, height =4, dpi = 1000)

bargraph04 <- finalSeurat %>%
  group_by(Timepoint, adj.ProtoOmi, res04Cluster) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))

ggplot(bargraph04, aes(fill = res04Cluster, y = Proportion, x = Timepoint))+
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Resolution 0.4")+
  facet_grid(cols = vars(adj.ProtoOmi))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(here::here("04_Analysis", "plots", "05_clustering", "Resolution04_ClustersOverTime.png"), width = 5, height =4, dpi = 1000)

DimPlot(seurat03res, split.by = "Timepoint", reduction = "harmony.wnn.umap")
DimPlot(seurat04res, split.by = "Timepoint", reduction = "harmony.wnn.umap")


#separate out by uninfected and infected
bargraph04 <- finalSeurat %>%
  group_by(Infection, Timepoint, adj.ProtoOmi, res04Cluster) %>%
  summarize(n = n()) %>%
  mutate(Proportion = n / sum(n),
         Timepoint = factor(Timepoint, levels = c("Day 0", "Day 15", "Day 90", "Day 180")))

ggplot(bargraph04, aes(fill = res04Cluster, y = Proportion, x = Timepoint))+
  geom_bar(position = "stack", stat = "identity")+
  ggtitle("Resolution 0.4")+
  facet_grid(cols = vars(adj.ProtoOmi),, rows = vars(Infection))+
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1, angle = 45))
ggsave(here::here("04_Analysis", "plots", "05_clustering", "Resolution04_ClustersOverTime_infection.png"), width = 5, height =7, dpi = 1000)
