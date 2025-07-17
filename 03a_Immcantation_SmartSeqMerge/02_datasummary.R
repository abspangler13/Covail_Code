#load dependencies
library(dplyr)
library(ggVennDiagram)
library(ggplot2)

#load in the data
df <- read.csv(here::here("03a_Immcantation_SmartSeqMerge", "MergedCiteSeqAndSmartSeqSequences.csv")) %>% mutate(clone_id = paste0(donor_id, "_", clone_id))

########what proportion of sequences in our CITESeq data have a clonemate from our Smartseq data - donut plot
stats <- df %>%
          group_by(clone_id, Dataset) %>%
          summarize(n = n()) %>%
          mutate(Overlap = length(unique(Dataset)) == 2,
                 sum = sum(n))

#make donut
donuts <- df %>%
            filter(Dataset == "CITESeq") %>%
            mutate(donut = case_when(clone_id %in% stats$clone_id[stats$Overlap == TRUE] ~"Overlapping",
                                     clone_id %in% stats$clone_id[stats$sum > 1] ~ "Non-overlapping",
                                     TRUE ~ "Singlet")) %>%
            mutate(clone_id = case_when(donut == "Singlet" ~ "Singlet",
                                        TRUE ~ clone_id)) %>%
            
            group_by(donut) %>%
            summarize(n = n()) %>%
            ungroup() %>%
            mutate(Proportion = n / sum(n)) %>%
            arrange(donut) %>%
            mutate(ymax = cumsum(Proportion),
                   ymin = c(0, head(ymax, n=-1)))
            
#plot donut!!!!
ggplot(donuts)+
  geom_rect(color= "black", linewidth=0, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=donut))+
  coord_polar(theta="y")+
  xlim(c(2,4))+
  scale_fill_manual(values = c( "Overlapping" = "firebrick", "Non-overlapping" = "pink", "Singlet" = "white"))+
  ggtitle("CITESeq Overlapping Clonal Groups")+
  theme_void()+
  theme(legend.title = element_blank())


#########make a venn diagram showing the overlap between both datasets in clonal families
x <- list(citeseq = unique(df$clone_id[df$Dataset == "CITESeq"]), smartseq = unique(df$clone_id[df$Dataset == "SmartSeq"]))
ggVennDiagram(x, label_size = 4, relative_width = 10, category.names = c("CITESeq", "SmartSeq"))+
  scale_x_continuous(expand = expansion(mult = .2))
ggsave(here::here("03a_Immcantation_SmartSeqMerge", "VennDiagramOverlap.png"), width = 7, height = 4)
