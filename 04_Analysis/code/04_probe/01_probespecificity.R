library(Seurat)
library(here)
library(dplyr)
library(tidyverse)
library(readr)
library(stringr)
library(ggplot2)
library(writexl)
library(readxl)
library(RColorBrewer)
library(tidyseurat)
library(devtools)

vdjObj <- readRDS(file = here::here("04_Analysis","data_objects","03_vdj","VDJPlusDSBNSeuratObj_NonVDJRemoved_CloneIDSubjPasted_final_all.rds"))

#write table describing counts per subject and timepoint
Cell.Subject <- vdjObj %>% group_by (Timepoint, Subject) %>% 
  dplyr::summarise(Freq = n()) %>% 
  pivot_wider(names_from = Timepoint, values_from = Freq)

write.csv(Cell.Subject,file = here::here("04_Analysis","data_objects","04_probe","Cells_per_timepoint_all.csv"))

#create scatterplots of all probes
probe <- str_replace_all(rownames(vdjObj@assays$Probes), "-", "_")
probe.combos <- combn(unique(probe),2)

pdf(file = here::here("04_Analysis","plots","04_probe","Probe_scatterplots_all.pdf"))
for(i in 1:ncol(probe.combos)){
  my.plot <- vdjObj %>% 
    join_features(features = c("Proto-RBD-PE", "BA1-RBD-PE", "XBB-RBD-no-fluor", "Beta-RBD-no-fluor"), assay="Probes") %>%
    pivot_wider(names_from=.feature,values_from=.abundance_Probes)
  
  colnames(my.plot) <- str_replace_all(colnames(my.plot), "-", "_")
  
  my.plot <- my.plot %>%
    ggplot((aes_string(x = probe.combos[1,i], y = probe.combos[2,i]))) + 
    geom_point(aes(color = Timepoint), size = 0.2, ratio = 1) + 
    theme_classic() +
    scale_x_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) + 
    scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 16, 32, 64, 128, 256, 5000)) +
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 10), 
          axis.title.y = element_text(size = 10), 
          axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8), 
          strip.text.x = element_text(size = 8)) + 
    facet_wrap(.~Timepoint) + 
    theme(aspect.ratio=1)
  plot(my.plot)
}
dev.off()

#here we will set thresholds for a probe positive label
#probe
# [1] "Proto_RBD_PE"      "BA1_RBD_PE"        "XBB_RBD_no_fluor"  "Beta_RBD_no_fluor"
threshold <- c(16, 18, 24, 16) # corresponds to indices above

#create jitterplots to determine thresholds for positivity
pdf(file = here::here("04_Analysis","plots","04_probe","JitterplotsForThresholdPositivity_ByOrigin.pdf"))
for(i in 1:length(probe)){
  my.plot <- vdjObj %>% 
    tidyseurat::join_features(features = c("Proto-RBD-PE", "BA1-RBD-PE", "XBB-RBD-no-fluor", "Beta-RBD-no-fluor"), assay="Probes") %>%
    pivot_wider(names_from=.feature,values_from=.abundance_Probes)
  
  colnames(my.plot) <- str_replace_all(colnames(my.plot), "-", "_")
  
  s <- my.plot %>%
    ggplot(aes(x=orig.ident, y=my.plot[[length(my.plot) - (length(probe) - i)]]))+
    geom_jitter(aes(color = orig.ident), size = 0.3, show.legend = FALSE)+
    geom_hline(yintercept=threshold[i])+
    labs(x="Sequencing Run", y=colnames(my.plot)[length(my.plot) - (length(probe) - i)], title=paste0(colnames(my.plot)[length(my.plot) - (length(probe) - i)], " Signal"))+
    scale_y_continuous(trans = "log2", limits = c(1, 5000), breaks = c(4, 8, 12, 16, 24, 32, 50, 128, 256, 5000)) +
    theme_classic()+
    theme(
      legend.position = 0, 
      axis.title.x = element_text(size = 10), 
      axis.title.y = element_text(size = 10), 
      axis.text.y = element_text(size = 8), 
      strip.text.x = element_text(size = 8),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=7))
  plot(s)
}
dev.off()
rm(my.plot)

#add a boolean value that describes probe positivity
probesColumns <- unlist(lapply(str_replace_all(probe,"_","-"), function(x) paste0(x, "_Positive")))

vdjDF <- vdjObj %>% 
  tidyseurat::join_features(features = c("Proto-RBD-PE", "BA1-RBD-PE", "XBB-RBD-no-fluor", "Beta-RBD-no-fluor"), assay="Probes") %>%
  pivot_wider(names_from=.feature,values_from=.abundance_Probes)

colnames(vdjDF) <- str_replace_all(colnames(vdjDF), "-", "_")

for(i in 1:length(probesColumns)){
  #define a positive value column name for each probe
  vdjDF[probesColumns[i]] <- vdjDF[probe[i]] >= threshold[(i)]
}

#set a temporary probe-positive pop to make some graphs
vdjDF$ProtoBA1 <- case_when(vdjDF$`Proto-RBD-PE_Positive` & vdjDF$`BA1-RBD-PE_Positive` ~ "Proto+BA1+",
                                !vdjDF$`Proto-RBD-PE_Positive` & vdjDF$`BA1-RBD-PE_Positive` ~ "Proto-BA1+",
                                vdjDF$`Proto-RBD-PE_Positive` & !vdjDF$`BA1-RBD-PE_Positive` ~ "Proto+BA1-",
                                TRUE ~ "Proto-BA1-")

pdf(file = here::here("04_Analysis", "plots", "04_probe", "ThresholdGatingOnly_ProtoBA1PopulationProportions_BeforeMergingBA1AndXBB.pdf"))
ggplot(vdjDF[vdjDF$ProtoBA1 != "Proto-BA1-",], aes(fill=ProtoBA1, y=1, x=Booster, color=ProtoBA1))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/BA1 Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster - Before Combining Omicron, Before gateR")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#combining XBB and BA1 for comparison after gateR
vdjDF$ProtoOmi <- case_when(vdjDF$`Proto-RBD-PE_Positive` & (vdjDF$`BA1-RBD-PE_Positive` | vdjDF$`XBB-RBD-no-fluor_Positive`) ~ "Proto+Omi+",
                                !vdjDF$`Proto-RBD-PE_Positive` & (vdjDF$`BA1-RBD-PE_Positive` | vdjDF$`XBB-RBD-no-fluor_Positive`) ~ "Proto-Omi+",
                                vdjDF$`Proto-RBD-PE_Positive` & !(vdjDF$`BA1-RBD-PE_Positive` | vdjDF$`XBB-RBD-no-fluor_Positive`) ~ "Proto+Omi-",
                                TRUE ~ "Proto-Omi-")

pdf(file = here::here("04_Analysis", "plots", "04_probe", "ThresholdGatingOnly_ProtoOmicronPopulationProportions_AfterMergingBA1AndXBB.pdf"))
ggplot(vdjDF[vdjDF$ProtoOmi != "Proto-Omi-",], aes(fill=ProtoOmi, y=1, x=Booster, color=ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/Omicron Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster - After Combining Omicron, Before gateR")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#make a scatterplot to show where the populations fall
pdf(file = here::here("04_Analysis", "plots", "04_probe", "ThresholdGatingOnly_ProtoBA1_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(vdjDF[vdjDF$ProtoOmi != "Proto-Omi-",], aes(y=BA1_RBD_PE, x=Proto_RBD_PE, color=ProtoOmi, fill=ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("BA1 RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs BA1 Colored By Population Labelling - XBB and BA1 Merged")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "04_probe", "ThresholdGatingOnly_ProtoXBB_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(vdjDF[vdjDF$ProtoOmi != "Proto-Omi-",], aes(y=XBB_RBD_no_fluor, x=Proto_RBD_PE, color=ProtoOmi, fill=ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("XBB RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs XBB Colored By Population Labelling - XBB and BA1 Merged")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()


#write a csv file that can be used downstream for gateR correction
write_csv(vdjDF[,c("X", "CELL","Booster", "Subject","Infection",colnames(vdjDF)[c(70:length(vdjDF))])], here::here("04_Analysis","data_objects", "04_probe", "thresholdgating_forcorrectionwithgateR.csv"))

#I will be putting clonal corrections in a separate R script so that we can correct clonal specificty
#calls separately- this allows us to flexibly choose which path we want

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
