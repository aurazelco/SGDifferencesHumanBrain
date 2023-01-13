main_Velm_3rd <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Velmeshev_2022_3rd_trimester/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/plot_high_variable_genes_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_3rd_trimester/cell_info_Velmeshev_2022_3rd_trimester.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)
#cell_info$ct <- str_replace_all(cell_info$ct, " ", "_")

ct_order <- c(
  "Dorsal progenitors",
  "Ventral progenitors",
  "Excitatory neurons",    
  "Interneurons",                      
  "OPCs",  
  "Oligodendrocytes",
  "Astrocytes",
  "Microglia",   
  "Vascular cells",  
  "Unknown"
)


pval_thresh <- 0.05
FC_thresh <- 1.2
runs <- c( "_1", "_2", "_3")

########## Plot cell types back to umap/tsne
  
Velm_3rd <- list()
for (v in runs) {
  Velm_3rd_input_seurat <- SCENICInput(main_Velm_3rd, F, v)
  Velm_3rd_seurat_list <- list()
  for (Velm_3rd_id in names(Velm_3rd_input_seurat)) {
    colnames(Velm_3rd_input_seurat[[Velm_3rd_id]]) <- str_replace_all(colnames(Velm_3rd_input_seurat[[Velm_3rd_id]]), "\\.", "-")
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(Velm_3rd_input_seurat[[Velm_3rd_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    Velm_3rd_seurat <- SCENICSeuratPlots(Velm_3rd_input_seurat, metadata_id, Velm_3rd_id)
    Velm_3rd_seurat_list <- append(Velm_3rd_seurat_list, list(Velm_3rd_seurat))
  }
  names(Velm_3rd_seurat_list) <- names(Velm_3rd_input_seurat)
  Velm_3rd <- append(Velm_3rd, list(Velm_3rd_seurat_list))
}
names(Velm_3rd) <- runs


k_clusters <- list("_1" = c(12, 12),
                   "_2" = c(12, 12),
                   "_3" = c(12, 12))

Velm_3rd_final <- list()
for (v in runs) { 
  Velm_3rd_seurat_list <- Velm_3rd[[v]]
  names(k_clusters[[v]]) <- names(Velm_3rd_seurat_list)
  for (Velm_3rd_id in names(Velm_3rd_seurat_list)) {
    #Velm_3rd_final[[v]][[Velm_3rd_id]] <- SCENICClustering(main_Velm_3rd, F, Velm_3rd_seurat_list[[Velm_3rd_id]], k_clusters[[v]][[Velm_3rd_id]], ct_order)
    # if also need to plot the UMAPs
    Velm_3rd_final[[v]][[Velm_3rd_id]] <- SCENICClustering(main_Velm_3rd, F, Velm_3rd_seurat_list[[Velm_3rd_id]], k_clusters[[v]][[Velm_3rd_id]], ct_order, plot_flag = "yes", Velm_3rd_id)
    SCENICMarkers(main_Velm_3rd, F, Velm_3rd[[v]][[Velm_3rd_id]], Velm_3rd_id)
  }
}

rm(Velm_3rd, Velm_3rd_id, Velm_3rd_input_seurat, Velm_3rd_seurat, Velm_3rd_seurat_list, k_clusters, metadata_id)

saveRDS(Velm_3rd_final, paste0(main_Velm_3rd, "seurat_files.rds"))

Velm_3rd_final <- readRDS(paste0(main_Velm_3rd, "seurat_files.rds"))

### Calculates Markers in each SeuratObject

Velm_3rd_markers <- SCENICInputMarkers(main_Velm_3rd, F, pval_thresh, FC_thresh)
Velm_3rd_10 <- SCENICtop10genes(Velm_3rd_markers, F)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Velm_3rd, F, Velm_3rd_final, Velm_3rd_10, ct_order)
#HmpSCENICAll(main_Velm_3rd, F, Velm_3rd_final, Velm_3rd_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

Velm_3rd_scenic <- ImportSCENICresults(main_Velm_3rd, F, "1_GRN", proj_order = "no")
#SCENICTfTg(main_Velm_3rd, F, Velm_3rd_scenic, Velm_3rd_final, ct_order)
SCENICTfTg(main_Velm_3rd, F, Velm_3rd_scenic, Velm_3rd_final, ct_order, 100)

Velm_3rd_tf_list <- SCENICExtractGRN(Velm_3rd_scenic, F, "TF", 100)
ExtractDiffGRN(main_Velm_3rd, F, Velm_3rd_tf_list, "TF")
SCENICPlotGRN(main_Velm_3rd, F, Velm_3rd_tf_list, "TF")

Velm_3rd_tg_list <- SCENICExtractGRN(Velm_3rd_scenic, F, "target", 50)
ExtractDiffGRN(main_Velm_3rd, F, Velm_3rd_tg_list, "Target")
SCENICPlotGRN(main_Velm_3rd, F, Velm_3rd_tg_list, "Target")

#####  TFs and Targets expression in original SeuratObject

Velm_3rd_trim <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_3rd_trimester.rds")
Velm_3rd_trim@meta.data$ct_sex <- paste(Velm_3rd_trim@meta.data$cluster_final, Velm_3rd_trim@meta.data$sex, sep="_")

Velm_3rd_trim_tf <- read.csv(paste0(main_Velm_3rd, "5_outputs/different_TF_between_sexes.csv"))
RidgeTFTG(main_Velm_3rd, Velm_3rd_trim, Velm_3rd_trim_tf$gene_id, "ct_sex", "TF")
RidgeTFTG(main_Velm_3rd, Velm_3rd_trim, Velm_3rd_trim_tf$gene_id, "sex", "TF")

Velm_3rd_trim_tg <- read.csv(paste0(main_Velm_3rd, "5_outputs/different_Target_between_sexes.csv"))
RidgeTFTG(main_Velm_3rd, Velm_3rd_trim, Velm_3rd_trim_tg$gene_id, "ct_sex", "Target")
RidgeTFTG(main_Velm_3rd, Velm_3rd_trim, Velm_3rd_trim_tg$gene_id, "sex", "Target")


#####  Regulons

# Velm_3rd
Velm_3rd_auc <- ImportSCENICresults(main_Velm_3rd, F, "3_AUCell", proj_order = "yes")
Velm_3rd_reg_list <- SCENICExtractRegulons(Velm_3rd_auc, F)
SCENICPlotRegulons(main_Velm_3rd, F, Velm_3rd_reg_list)

########## Number of TF-TG pairs between F and M of same project
Velm_3rd_overlapTFTG <- SCENICOverlapTfTg(Velm_3rd_scenic, F, analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_3rd, F, Velm_3rd_overlapTFTG)
Velm_3rd_overlapTFTG <- SCENICOverlapTfTg(Velm_3rd_scenic, F, threshold = 10000,  analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_3rd, F, Velm_3rd_overlapTFTG, threshold = 10000)


#####  Gene Variability

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_3rd_trimester/cell_info_Velmeshev_2022_3rd_trimester.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)

cell_info[which(cell_info$sex=="Female"), "sex"] <- "F"
cell_info[which(cell_info$sex=="Male"), "sex"] <- "M"

top2000 <- readRDS(paste0(main_Velm_3rd, "top_2000_SD_expr_matrix_Velmeshev_2022_3rd_trimester.rds"))

SexSD(main_Velm_3rd, cell_info, top2000)

