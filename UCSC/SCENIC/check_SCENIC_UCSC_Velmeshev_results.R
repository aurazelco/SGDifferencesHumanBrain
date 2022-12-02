main_Velm_2nd <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Velmeshev_2022_2nd_trimester/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/SCENIC/check_SCENIC_results_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2nd_trimester/cell_info_Velmeshev_2022_2nd_trimester.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)
#cell_info$ct <- str_replace_all(cell_info$ct, " ", "_")


ct_order <- c(
  "Dorsal progenitors",
  "Ventral progenitors",
  "Excitatory neurons",    
  "Interneurons",                      
  "OPCs",  
  "Astrocytes",
  "Microglia",   
  "Vascular cells",  
  "Unknown"
)

pval_thresh <- 0.05
FC_thresh <- 1.2
runs <- c( "_1", "_2", "_3")

########## Plot cell types back to umap/tsne

Velm_2nd <- list()
for (v in runs) {
  Velm_2nd_input_seurat <- SCENICInputSeurat(main_Velm_2nd, F, v)
  Velm_2nd_seurat_list <- list()
  for (Velm_2nd_id in names(Velm_2nd_input_seurat)) {
    colnames(Velm_2nd_input_seurat[[Velm_2nd_id]]) <- str_replace_all(colnames(Velm_2nd_input_seurat[[Velm_2nd_id]]), "\\.", "-")
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(Velm_2nd_input_seurat[[Velm_2nd_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    Velm_2nd_seurat <- SCENICSeuratPlots(Velm_2nd_input_seurat, metadata_id, Velm_2nd_id)
    Velm_2nd_seurat_list <- append(Velm_2nd_seurat_list, list(Velm_2nd_seurat))
  }
  names(Velm_2nd_seurat_list) <- names(Velm_2nd_input_seurat)
  Velm_2nd <- append(Velm_2nd, list(Velm_2nd_seurat_list))
}
names(Velm_2nd) <- runs


k_clusters <- list("_1" = c(12, 11),
                   "_2" = c(12, 11),
                   "_3" = c(13, 10))

Velm_2nd_final <- list()
for (v in runs) { 
  Velm_2nd_seurat_list <- Velm_2nd[[v]]
  names(k_clusters[[v]]) <- names(Velm_2nd_seurat_list)
  for (Velm_2nd_id in names(Velm_2nd_seurat_list)) {
    #Velm_2nd_final[[v]][[Velm_2nd_id]] <- SCENICClustering(main_Velm_2nd, F, Velm_2nd_seurat_list[[Velm_2nd_id]], k_clusters[[v]][[Velm_2nd_id]], ct_order)
    # if also need to plot the UMAPs
    Velm_2nd_final[[v]][[Velm_2nd_id]] <- SCENICClustering(main_Velm_2nd, F, Velm_2nd_seurat_list[[Velm_2nd_id]], k_clusters[[v]][[Velm_2nd_id]], ct_order, plot_flag = "yes", Velm_2nd_id)
    SCENICMarkers(main_Velm_2nd, F, Velm_2nd[[v]][[Velm_2nd_id]], Velm_2nd_id)
  }
}

rm(Velm_2nd, Velm_2nd_id, Velm_2nd_input_seurat, Velm_2nd_seurat, Velm_2nd_seurat_list, k_clusters, metadata_id)

saveRDS(Velm_2nd_final, paste0(main_Velm_2nd, "seurat_files.rds"))

Velm_2nd_final <- readRDS(paste0(main_Velm_2nd, "seurat_files.rds"))

### Calculates Markers in each SeuratObject

Velm_2nd_markers <- SCENICInputMarkers(main_Velm_2nd, F, pval_thresh, FC_thresh)
Velm_2nd_10 <- SCENICtop10genes(Velm_2nd_markers, F)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Velm_2nd, F, Velm_2nd_final, Velm_2nd_10, ct_order)
HmpSCENICAll(main_Velm_2nd, F, Velm_2nd_final, Velm_2nd_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

Velm_2nd_scenic <- SCENICresultsSeurat(main_Velm_2nd, F, "1_GRN", proj_order = "no")
SCENICTfTg(main_Velm_2nd, F, Velm_2nd_scenic, Velm_2nd_final, ct_order)
SCENICTfTg(main_Velm_2nd, F, Velm_2nd_scenic, Velm_2nd_final, ct_order, 100)

#####  Regulons

# Velm_2nd
Velm_2nd_auc <- SCENICresultsSeurat(main_Velm_2nd, F, "3_AUCell", proj_order = "yes")
Velm_2nd_reg_list <- SCENICExtractRegulons(Velm_2nd_auc, F)
SCENICPlotRegulons(main_Velm_2nd, F, Velm_2nd_reg_list)

########## Number of TF-TG pairs between F and M of same project
Velm_2nd_overlapTFTG <- SCENICOverlapTfTg(Velm_2nd_scenic, F, "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_2nd, F, Velm_2nd_overlapTFTG)

