library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)

disco_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/"
main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs/outputs/"

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))
Idents(disco_filt) <- "proj_sex_disease_ct"

final_groups <- read.csv(paste0(main, "final_filt.csv"))
final_groups[, 1] <- NULL
final_groups$idents <- paste(final_groups$proj, final_groups$sex, final_groups$disease, final_groups$ct, sep="_")
final_groups$name_subfolders <- str_replace_all(final_groups$ct, "/", "_")

col_factors <- c("proj", 
                 "sex",
                 "disease",
                 "ct",
                 "og",
                 "idents"
)

final_groups[col_factors] <- lapply(final_groups[col_factors], as.factor) 

sexes <- levels(final_groups$sex)

for (disease_type in levels(final_groups$disease)) {
  path <- paste0(main, disease_type, "01A_DEGs")
  dir.create(path, showWarnings = FALSE, recursive = T)
  final_filt_disease <- subset(final_groups, disease==disease_type)
  for (ct_type in unique(final_filt_disease$ct)) {
    if ((disease_type!="Multiple Sclerosis") & (length(unique(final_filt_disease[which(final_filt_disease$ct==ct_type), "proj"])) > 1)) {
      for (id in unique(final_filt_disease$proj)) {
        id1 <- paste(id, sexes[1], disease_type, ct_type, sep="_")
        id2 <- paste(id, sexes[2], disease_type, ct_type, sep="_")
        if ((id1 %in% levels(final_groups$idents)) & (id2 %in% levels(final_groups$idents))) {
          path_ct <- paste0(path, "/", unique(final_filt_disease[which(final_filt_disease$ct==ct_type), "name_subfolders"]))
          dir.create(path_ct, showWarnings = FALSE)
          deg1 <- FindMarkers(disco_filt, 
                                ident.1 = id1, 
                                ident.2 = id2,
                                logfc.threshold = 0.25,
                                min.pct = 0.1,
                                only.pos = TRUE)
          write.csv(deg1, paste0(path_ct, "/", id, "_", sexes[1], ".csv"))
          deg2 <- FindMarkers(disco_filt, 
                                ident.1 = id2, 
                                ident.2 = id1,
                                logfc.threshold = 0.25,
                                min.pct = 0.1,
                                only.pos = TRUE)
          write.csv(deg2, paste0(path_ct, "/", id, "_", sexes[2],".csv"))
        } 
      }
    } else if (disease_type=="Multiple Sclerosis") {
      for (id in unique(final_filt_disease$proj)) {
        id1 <- paste(id, sexes[1], disease_type, ct_type, sep="_")
        id2 <- paste(id, sexes[2], disease_type, ct_type, sep="_")
        if ((id1 %in% levels(final_groups$idents)) & (id2 %in% levels(final_groups$idents))) {
          path_ct <- paste0(path, "/", unique(final_filt_disease[which(final_filt_disease$ct==ct_type), "name_subfolders"]))
          dir.create(path_ct, showWarnings = FALSE)
          deg1 <- FindMarkers(disco_filt, 
                                ident.1 = id1, 
                                ident.2 = id2,
                                logfc.threshold = 0.25,
                                min.pct = 0.1,
                                only.pos = TRUE)
          write.csv(deg1, paste0(path_ct, "/", id, "_", sexes[1], ".csv"))
          deg2 <- FindMarkers(disco_filt, 
                                ident.1 = id2, 
                                ident.2 = id1,
                                logfc.threshold = 0.25,
                                min.pct = 0.1,
                                only.pos = TRUE)
          write.csv(deg2, paste0(path_ct, "/", id, "_", sexes[2],".csv"))
        } 
      }
    } else {
      print(paste(disease_type, ct_type, " has only one project with > 100 cells when there should be at least 2"))
    }
  }
}


write.csv(final_groups, paste0(main, "final_filt.csv"),row.names = F)
