library(Seurat)
library(tidyseurat)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(here)
library(stringr)
library(devtools)

#load gateR-corrected csv file
gateRData <- read.csv(here::here("04_Analysis", "data_objects", "04_probe", "gateRCorrected_Specificity.csv"))
#colnames(gateRData) <. str_replace_all(colnames(gateRData), "/.", "/.")

gateRData$BA1.RBD.PE_Positive <- gateRData$BA1_RBD_PE >= 4 #based on later functional data, 4 is a more appropriate threshold

#define the populations post.gateR correction
gateRData$Population <- case_when((gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/Beta+/BA1+/XBB+",
                              (gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/Beta+/BA1+",
                              (gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/BA1+",
                              (gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+",
                              (gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/Beta+",
                              (gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/XBB+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Beta+/BA1+/XBB+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Beta+/BA1+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "Beta+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Beta+/XBB+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "BA1+/XBB+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "XBB+",
                              !(gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & !gateRData$`XBB.RBD.no.fluor_Positive` ~ "BA1+",
                              (gateRData$`Proto.RBD.PE_Positive`) & (gateRData$`BA1.RBD.PE_Positive`) & !(gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/BA1+/XBB+",
                              (gateRData$`Proto.RBD.PE_Positive`) & !(gateRData$`BA1.RBD.PE_Positive`) & (gateRData$`Beta.RBD.no.fluor_Positive`) & gateRData$`XBB.RBD.no.fluor_Positive` ~ "Proto+/Beta+/XBB+",
                              TRUE ~ "RBD-")

gateRData$ProtoBA1 <- case_when(gateRData$Proto.RBD.PE_Positive & gateRData$BA1.RBD.PE_Positive ~ "Proto+BA1+",
                                !gateRData$Proto.RBD.PE_Positive & gateRData$BA1.RBD.PE_Positive ~ "Proto-BA1+",
                                gateRData$Proto.RBD.PE_Positive & !gateRData$BA1.RBD.PE_Positive ~ "Proto+BA1-",
                                TRUE ~ "Proto-BA1-")

pdf(file = here::here("04_Analysis", "plots", "04_probe", "gateRCorrected_ProtoBA1PopulationProportions_BeforeMergingBA1AndXBB.pdf"))
ggplot(gateRData[gateRData$ProtoBA1 != "Proto-BA1-",], aes(fill=ProtoBA1, y=1, x=Booster, color=ProtoBA1))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/BA1 Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster - Before Combining Omicron")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#combine the m'crons
gateRData$ProtoOmi <- case_when(gateRData$Proto.RBD.PE_Positive & (gateRData$BA1.RBD.PE_Positive | gateRData$XBB.RBD.no.fluor_Positive) ~ "Proto+Omi+",
                                                      !gateRData$Proto.RBD.PE_Positive & (gateRData$BA1.RBD.PE_Positive | gateRData$XBB.RBD.no.fluor_Positive) ~ "Proto-Omi+",
                                                      gateRData$Proto.RBD.PE_Positive & !(gateRData$BA1.RBD.PE_Positive | gateRData$XBB.RBD.no.fluor_Positive) ~ "Proto+Omi-",
                                                      TRUE ~ "Proto-Omi-")

pdf(file = here::here("04_Analysis", "plots", "04_probe", "gateRCorrected_ProtoOmicronPopulationProportions_AfterMergingBA1AndXBB.pdf"))
ggplot(gateRData[gateRData$ProtoOmi != "Proto-Omi-",], aes(fill=ProtoOmi, y=1, x=Booster, color=ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/Omicron Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions per Booster - After Combining Omicron")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#let's make a plot that shows where the populations lie
pdf(file = here::here("04_Analysis", "plots", "04_probe", "gateRCorrectedData_ProtoBA1_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(gateRData[gateRData$ProtoOmi != "Proto-Omi-",], aes(y=BA1_RBD_PE, x=Proto_RBD_PE, color=ProtoOmi, fill=ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("BA1 RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs BA1 Colored By Population Labelling - gateR Corrected and XBB and BA1 Merged")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "04_probe", "gateRCorrectedData_ProtoXBB_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(gateRData[gateRData$ProtoOmi != "Proto-Omi-",], aes(y=XBB_RBD_no_fluor, x=Proto_RBD_PE, color=ProtoOmi, fill=ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("XBB RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs XBB Colored By Population Labelling - gateR Corrected and XBB and BA1 Merged")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

######merge the new annotations with the seurat object metadata
vdjObj <- readRDS(file = here::here("04_Analysis","data_objects","03_vdj","VDJPlusDSBNSeuratObj_NonVDJRemoved_CloneIDSubjPasted_final_all.rds"))

counts <- vdjObj %>% 
  tidyseurat::join_features(features = c("Proto-RBD-PE", "BA1-RBD-PE", "XBB-RBD-no-fluor", "Beta-RBD-no-fluor"), assay="Probes") %>%
  pivot_wider(names_from=.feature,values_from=.abundance_Probes)

counts <- as.data.frame(counts)

rownames(counts) <- counts$CELL
rownames(gateRData) <- gateRData$CELL

#add in HTO signal data
#set old object
oldObj <- vdjObj@meta.data
vdjObj <- AddMetaData(vdjObj, as.data.frame(counts[,colnames(counts) %in% c("Proto-RBD-PE", "BA1-RBD-PE", "XBB-RBD-no-fluor", "Beta-RBD-no-fluor")]))
vdjObj <- vdjObj %>% filter(`BA1-RBD-PE` <= 2000)

#add in positivity booleans
vdjObj <- AddMetaData(vdjObj, gateRData[,colnames(gateRData) %in% c("Proto.RBD.PE_Positive","BA1.RBD.PE_Positive", "XBB.RBD.no.fluor_Positive", "Beta.RBD.no.fluor_Positive", "Population", "ProtoBA1", "ProtoOmi")])


#240119: Adding back in CD20 to the prot data brought in new cells not present before, and by gating they all look cross reactive, so I will label them as such here
#In the future, we might want to fix this more formally by gating again in gateR, but for now, it's not worth tinkering too much with the gating

#if there are any notes manually added via gateR, then add them in
if("notes" %in% colnames(gateRData)){
  vdjObj <- AddMetaData(vdjObj, gateRData[,colnames(gateRData) %in% c("notes")])
}

write.csv(vdjObj@meta.data, file=here::here("04_Analysis", "data_objects", "04_probe", "COVAILVDJAndCSOData_SpecificitiesLabelled_ClonesNotCorrected.csv"))

#########Set the clonal group's specificity so that it is the same for all of them
#start with general population labels
vdjObj$CloneSubjectIDTimepoint <- paste0(vdjObj$clone_subject_id, "_", vdjObj$Timepoint)
vdjDF <- vdjObj@meta.data

adj.spec.P <- as.data.frame(vdjDF %>%
                              group_by(clone_subject_id,Population) %>% 
                              dplyr::summarise(Freq = n()) %>% 
                              pivot_wider(names_from = Population, values_from = Freq) %>% 
                              replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$clone_subject_id #this block takes the summary table and chooses the max label as the clonal specificity
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P$adj.Population<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("adj.Population","clone_subject_id")]
vdjDF <- vdjDF %>% left_join(adj.spec.P,by="clone_subject_id")

#combine with protoomi now
adj.spec.P <- as.data.frame(vdjDF %>%
                              group_by(clone_subject_id,ProtoOmi) %>% 
                              dplyr::summarise(Freq = n()) %>% 
                              pivot_wider(names_from = ProtoOmi, values_from = Freq) %>% 
                              replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$clone_subject_id
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P$adj.ProtoOmi<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$clone_subject_id <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("adj.ProtoOmi","clone_subject_id")]
vdjDF <- vdjDF %>% left_join(adj.spec.P,by="clone_subject_id")

#Make a second adjusted specificity label that is relative to a timepoint - this will capture potential changes in specificity over time
adj.spec.P <- as.data.frame(vdjDF %>%
                              group_by(CloneSubjectIDTimepoint,ProtoOmi) %>% 
                              dplyr::summarise(Freq = n()) %>% 
                              pivot_wider(names_from = ProtoOmi, values_from = Freq) %>% 
                              replace(is.na(.),0))

rownames(adj.spec.P) <- adj.spec.P$CloneSubjectIDTimepoint
adj.spec.P <- adj.spec.P[,-1]
adj.spec.P$TimeAdj.ProtoOmi<-colnames(adj.spec.P)[apply(adj.spec.P,1,which.max)]
adj.spec.P$CloneSubjectIDTimepoint <- rownames(adj.spec.P)
adj.spec.P <- adj.spec.P[,c("TimeAdj.ProtoOmi","CloneSubjectIDTimepoint")]
vdjDF <- vdjDF %>% left_join(adj.spec.P,by="CloneSubjectIDTimepoint")

#let's plot the adjusted value proportions
pdf(file = here::here("04_Analysis", "plots", "04_probe", "ClonallyAdjusted_gateRCorrected_ProtoOmicronPopulationProportions_AfterMergingBA1AndXBB.pdf"))
ggplot(vdjDF[vdjDF$adj.ProtoOmi != "Proto-Omi-",], aes(fill=adj.ProtoOmi, y=1, x=Booster, color=adj.ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/Omicron Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions - Omicrons Combined, Clones Corrected")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#let's see where they sit on the proto-omi plots
pdf(file = here::here("04_Analysis", "plots", "04_probe", "ClonallyAdjusted_gateRCorrectedData_ProtoBA1_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(vdjDF[vdjDF$adj.ProtoOmi != "Proto-Omi-",], aes(y=`BA1-RBD-PE`, x=`Proto-RBD-PE`, color=adj.ProtoOmi, fill=adj.ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("BA1 RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs BA1 Colored By Population Labelling - ClonallyCorrected")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "04_probe", "ClonallyAdjusted_gateRCorrectedData_ProtoXBB_BA1XBBMergedLabelledScatterplot.pdf"))
ggplot(vdjDF[vdjDF$adj.ProtoOmi != "Proto-Omi-",], aes(y=`XBB-RBD-no-fluor`, x=`Proto-RBD-PE`, color=adj.ProtoOmi, fill=adj.ProtoOmi))+
  geom_point()+
  scale_x_continuous(trans='log2')+
  scale_y_continuous(trans='log2')+
  ylab("XBB RBD Counts")+
  xlab("Prototype RBD Counts")+
  ggtitle("Proto vs XBB Colored By Population Labelling - ClonallyCorrected")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), axis.text.y=element_text(size=12))
dev.off()

#let's make one last comparison- compare overall proportions between boosters + infection status
vdjDF$BoostInfect <- paste0(vdjDF$Booster, "_", vdjDF$Infection)

pdf(file = here::here("04_Analysis", "plots", "04_probe", "AdjustedPopulations_ByBoosterAndInfectionStatus.pdf"))
ggplot(vdjDF[vdjDF$adj.ProtoOmi != "Proto-Omi-",], aes(fill=adj.ProtoOmi, y=1, x=BoostInfect, color=adj.ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/Omicron Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions - Adjusted Spec.")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90), axis.text.y=element_text(size=12))
dev.off()

pdf(file = here::here("04_Analysis", "plots", "04_probe", "TimeAdjustedPopulations_ByBoosterAndInfectionStatus.pdf"))
ggplot(vdjDF[vdjDF$TimeAdj.ProtoOmi != "Proto-Omi-",], aes(fill=TimeAdj.ProtoOmi, y=1, x=BoostInfect, color=TimeAdj.ProtoOmi))+
  geom_bar(position="fill",stat="identity")+
  ylab("Percentage of Booster Proto/Omicron Populations")+
  xlab("Booster")+
  ggtitle("Probe+ Proportions - Time-Adjusted Spec.")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90), axis.text.y=element_text(size=12))
dev.off()

#append metadata
rownames(vdjDF) <- vdjDF$CELL

vdjObj <- AddMetaData(vdjObj, vdjDF[,colnames(vdjDF) %in% c("adj.Population", "adj.ProtoOmi", "TimeAdj.ProtoOmi")])

#########write new seurat object and export metadata for Flavio to use
saveRDS(vdjObj, file=here::here("04_Analysis", "data_objects", "04_probe", "CoVSeuratObj_VDJCSOGEX_SpecificitiesLabelled_CloneCorrected.rds"))

write.csv(vdjObj@meta.data, file=here::here("04_Analysis", "data_objects", "04_probe", "COVAILVDJAndCSOData_SpecificitiesLabelled_ClonesCorrected.csv"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

