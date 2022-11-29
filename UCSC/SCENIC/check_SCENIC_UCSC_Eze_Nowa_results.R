main_Eze_Nowa <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Eze_Nowakowski_integrated_2nd_trimester/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Eze_Nowakowski_integrated_2nd_trimester/cell_info.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)
#cell_info$ct <- str_replace_all(cell_info$ct, " ", "_")

ct_order <- c(
  "Mesenchymal",
  "Neuroepithelial",
  "Neuronal",
  "Radial Glial",
  "Other"
)

pval_thresh <- 0.05
FC_thresh <- 1.2
runs <- c( "_1", "_2", "_3")

########## Plot cell types back to umap/tsne

eze_nowa <- list()
for (v in runs) {
  eze_nowa_input_seurat <- SCENICInputSeurat(main_Eze_Nowa, F, v)
  eze_nowa_seurat_list <- list()
  for (eze_nowa_id in names(eze_nowa_input_seurat)) {
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(eze_nowa_input_seurat[[eze_nowa_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    eze_nowa_seurat <- SCENICSeuratPlots(eze_nowa_input_seurat, metadata_id, eze_nowa_id)
    eze_nowa_seurat_list <- append(eze_nowa_seurat_list, list(eze_nowa_seurat))
  }
  names(eze_nowa_seurat_list) <- names(eze_nowa_input_seurat)
  eze_nowa <- append(eze_nowa, list(eze_nowa_seurat_list))
}
names(eze_nowa) <- runs

k_clusters <- list("_1" = c(12, 9),
                   "_2" = c(12, 9),
                   "_3" = c(12, 9))

eze_nowa_final <- list()
for (v in runs) { 
  eze_nowa_seurat_list <- eze_nowa[[v]]
  names(k_clusters[[v]]) <- names(eze_nowa_seurat_list)
  for (eze_nowa_id in names(eze_nowa_seurat_list)) {
    #eze_nowa_final[[v]][[eze_nowa_id]] <- SCENICClustering(main_Eze_Nowa, F, eze_nowa_seurat_list[[eze_nowa_id]], k_clusters[[v]][[eze_nowa_id]], ct_order)
    # if also need to plot the UMAPs
    eze_nowa_final[[v]][[eze_nowa_id]] <- SCENICClustering(main_Eze_Nowa, F, eze_nowa_seurat_list[[eze_nowa_id]], k_clusters[[v]][[eze_nowa_id]], ct_order, plot_flag = "yes", eze_nowa_id)
    SCENICMarkers(main_Eze_Nowa, F, eze_nowa[[v]][[eze_nowa_id]], eze_nowa_id)
  }
}

rm(eze_nowa, eze_nowa_id, eze_nowa_input_seurat, eze_nowa_seurat, eze_nowa_seurat_list, k_clusters, metadata_id)

saveRDS(eze_nowa_final, paste0(main_Eze_Nowa, "seurat_files.rds"))

eze_nowa_final <- readRDS(paste0(main_Eze_Nowa, "seurat_files.rds"))

### Calculates Markers in each SeuratObject

eze_nowa_markers <- SCENICInputMarkers(main_Eze_Nowa, F, pval_thresh, FC_thresh)
eze_nowa_10 <- SCENICtop10genes(eze_nowa_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Eze_Nowa, F, eze_nowa_final, eze_nowa_10, ct_order)
HmpSCENICAll(main_Eze_Nowa, F, eze_nowa_final, eze_nowa_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

eze_nowa_scenic <- SCENICresultsSeurat(main_Eze_Nowa, F, "1_GRN", proj_order = "no")
SCENICTfTg(main_Eze_Nowa, F, eze_nowa_scenic, eze_nowa_final, ct_order)
SCENICTfTg(main_Eze_Nowa, F, eze_nowa_scenic, eze_nowa_final, ct_order, 100)

#####  Regulons

# eze_nowa
eze_nowa_auc <- SCENICresultsSeurat(main_Eze_Nowa, F, "3_AUCell", proj_order = "yes")
eze_nowa_reg_list <- SCENICExtractRegulons(eze_nowa_auc, F)
SCENICPlotRegulons(main_Eze_Nowa, F, eze_nowa_reg_list)

########## Number of TF-TG pairs between F and M of same project
eze_nowa_overlapTFTG <- SCENICOverlapTfTg(eze_nowa_scenic, F)
SCENICPlotOverlapTfTg(main_Eze_Nowa, F, eze_nowa_overlapTFTG)

