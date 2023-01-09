main_Velm_2nd_year <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Velmeshev_2022_1_2_years/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/SCENIC/check_SCENIC_results_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_1_2_years/cell_info_Velmeshev_2022_1_2_years.csv")
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

Velm_2nd_year <- list()
for (v in runs) {
  Velm_2nd_year_input_seurat <- SCENICInputSeurat(main_Velm_2nd_year, F, v)
  Velm_2nd_year_seurat_list <- list()
  for (Velm_2nd_year_id in names(Velm_2nd_year_input_seurat)) {
    colnames(Velm_2nd_year_input_seurat[[Velm_2nd_year_id]]) <- str_replace_all(colnames(Velm_2nd_year_input_seurat[[Velm_2nd_year_id]]), "\\.", "-")
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(Velm_2nd_year_input_seurat[[Velm_2nd_year_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    Velm_2nd_year_seurat <- SCENICSeuratPlots(Velm_2nd_year_input_seurat, metadata_id, Velm_2nd_year_id)
    Velm_2nd_year_seurat_list <- append(Velm_2nd_year_seurat_list, list(Velm_2nd_year_seurat))
  }
  names(Velm_2nd_year_seurat_list) <- names(Velm_2nd_year_input_seurat)
  Velm_2nd_year <- append(Velm_2nd_year, list(Velm_2nd_year_seurat_list))
}
names(Velm_2nd_year) <- runs


k_clusters <- list("_1" = c(10, 12),
                   "_2" = c(11, 12),
                   "_3" = c(10, 12))

Velm_2nd_year_final <- list()
for (v in runs) { 
  Velm_2nd_year_seurat_list <- Velm_2nd_year[[v]]
  names(k_clusters[[v]]) <- names(Velm_2nd_year_seurat_list)
  for (Velm_2nd_year_id in names(Velm_2nd_year_seurat_list)) {
    #Velm_2nd_year_final[[v]][[Velm_2nd_year_id]] <- SCENICClustering(main_Velm_2nd_year, F, Velm_2nd_year_seurat_list[[Velm_2nd_year_id]], k_clusters[[v]][[Velm_2nd_year_id]], ct_order)
    # if also need to plot the UMAPs
    Velm_2nd_year_final[[v]][[Velm_2nd_year_id]] <- SCENICClustering(main_Velm_2nd_year, F, Velm_2nd_year_seurat_list[[Velm_2nd_year_id]], k_clusters[[v]][[Velm_2nd_year_id]], ct_order, plot_flag = "yes", Velm_2nd_year_id)
    SCENICMarkers(main_Velm_2nd_year, F, Velm_2nd_year[[v]][[Velm_2nd_year_id]], Velm_2nd_year_id)
  }
}

rm(Velm_2nd_year, Velm_2nd_year_id, Velm_2nd_year_input_seurat, Velm_2nd_year_seurat, Velm_2nd_year_seurat_list, k_clusters, metadata_id)

saveRDS(Velm_2nd_year_final, paste0(main_Velm_2nd_year, "seurat_files.rds"))

Velm_2nd_year_final <- readRDS(paste0(main_Velm_2nd_year, "seurat_files.rds"))

### Calculates Markers in each SeuratObject

Velm_2nd_year_markers <- SCENICInputMarkers(main_Velm_2nd_year, F, pval_thresh, FC_thresh)
Velm_2nd_year_10 <- SCENICtop10genes(Velm_2nd_year_markers, F)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Velm_2nd_year, F, Velm_2nd_year_final, Velm_2nd_year_10, ct_order)
#HmpSCENICAll(main_Velm_2nd_year, F, Velm_2nd_year_final, Velm_2nd_year_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

Velm_2nd_year_scenic <- SCENICresultsSeurat(main_Velm_2nd_year, F, "1_GRN", proj_order = "no")
#SCENICTfTg(main_Velm_2nd_year, F, Velm_2nd_year_scenic, Velm_2nd_year_final, ct_order)
SCENICTfTg(main_Velm_2nd_year, F, Velm_2nd_year_scenic, Velm_2nd_year_final, ct_order, 100)

Velm_2nd_year_tf_list <- SCENICExtractGRN(Velm_2nd_year_scenic, F, "TF", 100)
ExtractDiffGRN(main_Velm_2nd_year, F, Velm_2nd_year_tf_list, "TF")
SCENICPlotGRN(main_Velm_2nd_year, F, Velm_2nd_year_tf_list, "TF")

Velm_2nd_year_tg_list <- SCENICExtractGRN(Velm_2nd_year_scenic, F, "target", 50)
ExtractDiffGRN(main_Velm_2nd_year, F, Velm_2nd_year_tg_list, "Target")
SCENICPlotGRN(main_Velm_2nd_year, F, Velm_2nd_year_tg_list, "Target")


#####  TFs and Targets expression in original SeuratObject

Velm_2nd_year <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_1_2_years.rds")
Velm_2nd_year@meta.data$ct_sex <- paste(Velm_2nd_year@meta.data$cluster_final, Velm_2nd_year@meta.data$sex, sep="_")

Velm_2nd_year_tf <- read.csv(paste0(main_Velm_2nd_year, "5_outputs/different_TF_between_sexes.csv"))
RidgeTFTG(main_Velm_2nd_year, Velm_2nd_year, Velm_2nd_year_tf$gene_id, "ct_sex", "TF")

Velm_2nd_year_tg <- read.csv(paste0(main_Velm_2nd_year, "5_outputs/different_Target_between_sexes.csv"))
RidgeTFTG(main_Velm_2nd_year, Velm_2nd_year, Velm_2nd_year_tg$gene_id, "ct_sex", "Target")


#####  Regulons

# Velm_2nd_year
Velm_2nd_year_auc <- SCENICresultsSeurat(main_Velm_2nd_year, F, "3_AUCell", proj_order = "yes")
Velm_2nd_year_reg_list <- SCENICExtractRegulons(Velm_2nd_year_auc, F)
SCENICPlotRegulons(main_Velm_2nd_year, F, Velm_2nd_year_reg_list)

########## Number of TF-TG pairs between F and M of same project
Velm_2nd_year_overlapTFTG <- SCENICOverlapTfTg(Velm_2nd_year_scenic, F, "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_2nd_year, F, Velm_2nd_year_overlapTFTG)

#####  Gene Variability

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/SCENIC/plot_high_variable_genes_func.R")

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_1_2_years/cell_info_Velmeshev_2022_1_2_years.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)

cell_info[which(cell_info$sex=="Female"), "sex"] <- "F"
cell_info[which(cell_info$sex=="Male"), "sex"] <- "M"

top2000 <- readRDS(paste0(main_Velm_2nd_year, "top_2000_SD_expr_matrix_Velmeshev_2022_1_2_years.rds"))

SexSD(main_Velm_2nd_year, cell_info, top2000)

