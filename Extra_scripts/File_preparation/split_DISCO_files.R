library(tidyr)
library(stringr)

main_out <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/HumanBrainSexSingleCell/data/DEGs/"
tot_genes <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/tot_genes_ct.csv")

tot_genes$X <- NULL

tot_genes$disease <- str_replace_all(tot_genes$disease, "Normal", "Healthy")

tot_genes_ls <- split(tot_genes, f = tot_genes$disease)

lapply(1:length(tot_genes_ls), function(x) write.csv(tot_genes_ls[[x]], paste0(main_out, "tot_genes_", names(tot_genes_ls)[x], ".csv")))


cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/cell_info.csv")
cell_info$X <- NULL
cell_info$og_group <- str_replace_all(cell_info$og_group, "Normal", "Healthy")

cell_info <- separate(cell_info, og_group, into = c("proj", "sex", "disease", "ct"), sep = "_", remove = F)


cell_info$proj_dis <- paste(cell_info$proj, cell_info$disease, sep="_")
cell_info <- cell_info[, c(1,2,3,4,9)]
cell_info_ls <- split(cell_info, f = cell_info$proj_dis)
for (i in names(cell_info_ls)) {
  cell_info_ls[[i]] <- cell_info_ls[[i]][, c(1,2,3,4)]
}

lapply(1:length(cell_info_ls), function(x) write.csv(cell_info_ls[[x]], paste0(main_out, "cell_info_", names(cell_info_ls)[x], ".csv")))