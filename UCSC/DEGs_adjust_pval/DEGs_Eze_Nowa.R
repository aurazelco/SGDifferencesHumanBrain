# 01A_generate_DEGs.R

library(Seurat)

main <- "UCSC/DEGs_2nd_trimester/Eze_Nowakowski_integrated_2nd_trimester/outputs/"

trim_2nd <- readRDS("UCSC/Seurat_UCSC/Eze_Nowakowski_integrated_2nd_trimester.rds")
Idents(trim_2nd) <- "sex_ct"

min_cells <- 100

final_groups <- read.csv(paste0(main, "final_filt_", min_cells, ".csv"))
sexes <- unique(final_groups$sex)

path <- paste0(main, "01A_DEGs")
dir.create(path, showWarnings = FALSE, recursive = T)
for (ct_type in unique(final_groups$ct)) {
  id1 <- paste(sexes[1], ct_type, sep="_")
  id2 <- paste(sexes[2], ct_type, sep="_")
  if ((id1 %in% unique(final_groups$idents)) & (id2 %in% unique(final_groups$idents))) {
    path_ct <- paste0(path, "/", ct_type)
    dir.create(path_ct, showWarnings = FALSE)
    deg1 <- FindMarkers(trim_2nd, 
                              ident.1 = id1, 
                              ident.2 = id2,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
    write.csv(deg1, paste0(path_ct, "/", sexes[1], ".csv"))
    deg2 <- FindMarkers(trim_2nd, 
                              ident.1 = id2, 
                              ident.2 = id1,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
    write.csv(deg2, paste0(path_ct, "/", sexes[2],".csv"))
  }
}


