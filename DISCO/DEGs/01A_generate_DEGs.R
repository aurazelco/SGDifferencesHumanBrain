setwd("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/")
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)


disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")
disco_filt@meta.data$proj_sex_disease_ct <- lapply(disco_filt@meta.data$proj_sex_disease_ct, str_replace_all, pattern="/", replacement="_")
Idents(disco_filt) <- "proj_sex_disease_ct"


final_groups <- read.csv("final_filt.csv")
final_groups[, 1] <- NULL

# had to do this because / in ct names were causing issues with dir.create
final_groups$ct <- sapply(final_groups$ct, str_replace_all, pattern="/", replacement="_")
final_groups$idents <- sapply(final_groups$idents, str_replace_all, pattern="/", replacement="_")

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
  path <- paste0(getwd(), "/", disease_type)
  dir.create(path, showWarnings = FALSE)
  for (ct_type in levels(final_groups$ct)) {
      for (id in levels(final_groups$proj)) {
        id1 <- paste(id, sexes[1], disease_type, ct_type, sep="_")
        id2 <- paste(id, sexes[2], disease_type, ct_type, sep="_")
        if ((id1 %in% levels(final_groups$idents)) & (id2 %in% levels(final_groups$idents))) {
          path <- paste0(getwd(), "/", disease_type, "/", ct_type)
          dir.create(path, showWarnings = FALSE)
          deg1 <- FindMarkers(disco_filt, 
                             ident.1 = id1, 
                             ident.2 = id2,
                             logfc.threshold = 0.25,
                             min.pct = 0.1,
                             only.pos = TRUE)
          write.csv(deg1, paste0(path, "/", id, "_", sexes[1], ".csv"))
          deg2 <- FindMarkers(disco_filt, 
                              ident.1 = id2, 
                              ident.2 = id1,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
          write.csv(deg2, paste0(path, "/", id, "_", sexes[2],".csv"))
        }
      }
  }
}



