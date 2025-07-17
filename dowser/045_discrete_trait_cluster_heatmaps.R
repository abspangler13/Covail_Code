library(pheatmap)
library(dplyr)
library(tidyr)
library(textshape)
library(RColorBrewer)

dat <- read.csv(file = here::here("analysis","data_objects","14_dowser_redo","testSP_cluster_all.csv"), row.names=1)
# dat <- read.csv(file = here::here("analysis","data_objects","14_dowser_redo","testSP_cluster_all_connection.csv"), row.names=1)

dat <- dat %>% select(TO,FROM,PGT)
dat_wide <- dat %>% pivot_wider(names_from = TO, values_from = PGT)
dat_matrix <- dat_wide %>%
  column_to_rownames(loc = "FROM") %>%
  as.matrix()
dat_matrix[is.na(dat_matrix)] <- 1

##rename columns
colnames(dat_matrix)[colnames(dat_matrix) == "AM1a-6"] <- "Act_C4"
colnames(dat_matrix)[colnames(dat_matrix) == "AM2-0"] <- "Int_C3"
colnames(dat_matrix)[colnames(dat_matrix) == "AM3-acute-1"] <- "Atyp_C6"
colnames(dat_matrix)[colnames(dat_matrix) == "RM-2"] <- "Rest_C2"
colnames(dat_matrix)[colnames(dat_matrix) == "IgM-5"] <- "Unswitch_C1"
colnames(dat_matrix)[colnames(dat_matrix) == "AM3-chronic-3"] <- "Atyp_C7"
colnames(dat_matrix)[colnames(dat_matrix) == "AM1a-4"] <- "Act_C5"
## rename rows
rownames(dat_matrix)[rownames(dat_matrix) == "AM1a-4"] <- "Act_C5"
rownames(dat_matrix)[rownames(dat_matrix) == "AM3-chronic-3"] <- "Atyp_C7"
rownames(dat_matrix)[rownames(dat_matrix) == "RM-2"] <- "Rest_C2"
rownames(dat_matrix)[rownames(dat_matrix) == "IgM-5"] <- "Unswitch_C1"
rownames(dat_matrix)[rownames(dat_matrix) == "AM1a-6"] <- "Act_C4"
rownames(dat_matrix)[rownames(dat_matrix) == "AM2-0"] <- "Int_C3"
rownames(dat_matrix)[rownames(dat_matrix) == "AM3-acute-1"] <- "Atyp_C6"

#reorder columns
dat_matrix <- dat_matrix[c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7"),c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7")]
# dat_matrix <- dat_matrix[c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6"),c("Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7")]

# Create a heatmap
pdf(file = here::here("analysis", "plots", "14_dowser_redo", "heatmap_cluster_all.pdf"))

pheatmap(
  dat_matrix, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  display_numbers = TRUE,
  number_format = "%.2f", 
  number_color = "black", 
  fontsize_number = 20,
  color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100),
  breaks = c(
    seq(min(dat_matrix), 0.05, length.out = 50), 
    seq(0.051, max(dat_matrix), length.out = 50)
  ),
  fontsize_row = 15, 
  fontsize_col = 15, 
  legend = TRUE, 
  legend_labels = list(fontsize = 15)
)

dev.off()


## specificity exposure
dat <- read.csv(file = here::here("analysis","data_objects","14_dowser_redo","testSP_cluster_specificity_exposure.csv"))

# dat <- read.csv(file = here::here("analysis","data_objects","14_dowser_redo","testSP_cluster_specificity_exposure_connection.csv"))

dat_list <- split(dat, dat$specificity_exposure)

dat_matrices <- lapply(dat_list, function(sub_dat) {
    sub_dat <- sub_dat %>% select(TO, FROM, PGT)
    sub_dat_wide <- sub_dat %>% pivot_wider(names_from = TO, values_from = PGT)
    sub_dat_matrix <- sub_dat_wide %>%
        column_to_rownames(loc = "FROM") %>%
        as.matrix()
    sub_dat_matrix[is.na(sub_dat_matrix)] <- 1
    return(sub_dat_matrix)
})

# ## add in missing rows and columns
# dat_matrices$H2_Cross_Naive <- rbind(dat_matrices$H2_Cross_Naive, `IgM-5` = rep(1, ncol(dat_matrices$H2_Cross_Naive)))
# dat_matrices$H2_Cross_Naive <- cbind(dat_matrices$H2_Cross_Naive, `IgM-5` = rep(1, nrow(dat_matrices$H2_Cross_Naive)))

for (name in names(dat_matrices)) {
    #reorder rows and columns
    dat_matrices[[name]] <- dat_matrices[[name]][c("IgM-5","RM-2","AM2-0","AM1a-6","AM1a-4","AM3-acute-1","AM3-chronic-3"), c("IgM-5","RM-2","AM2-0","AM1a-6","AM1a-4","AM3-acute-1","AM3-chronic-3")]
    #rename rows and columns
    colnames(dat_matrices[[name]]) <- c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7")
    rownames(dat_matrices[[name]]) <- c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7")

    ##rename columns
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "AM1a-6"] <- "Act_C4"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "AM2-0"] <- "Int_C3"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "AM3-acute-1"] <- "Atyp_C6"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "RM-2"] <- "Rest_C2"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "IgM-5"] <- "Unswitch_C1"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "AM3-chronic-3"] <- "Atyp_C7"
    colnames(dat_matrices[[name]])[colnames(dat_matrices[[name]]) == "AM1a-4"] <- "Act_C5"
    ## rename rows
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "AM1a-4"] <- "Act_C5"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "AM3-chronic-3"] <- "Atyp_C7"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "RM-2"] <- "Rest_C2"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "IgM-5"] <- "Unswitch_C1"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "AM1a-6"] <- "Act_C4"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "AM2-0"] <- "Int_C3"
    rownames(dat_matrices[[name]])[rownames(dat_matrices[[name]]) == "AM3-acute-1"] <- "Atyp_C6"

#reorder columns
    dat_matrices[[name]] <- dat_matrices[[name]][c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7"),c("Unswitch_C1","Rest_C2","Int_C3","Act_C4","Act_C5","Atyp_C6","Atyp_C7")]
}
for (name in names(dat_matrices)) {
    pdf(file = here::here("analysis", "plots", "14_dowser_redo", paste0("heatmap_cluster_", name, ".pdf")))
    pheatmap(dat_matrices[[name]], cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, 
             number_format = "%.2f", number_color = "black", fontsize_number = 20,
             color = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100), 
             breaks = c(seq(min(dat_matrices[[name]]), 0.05, length.out = 50), seq(0.051, max(dat_matrices[[name]]), length.out = 50)),
             fontsize_row = 18, fontsize_col = 18, legend = TRUE, legend_labels = list(fontsize = 15))
    dev.off()
}