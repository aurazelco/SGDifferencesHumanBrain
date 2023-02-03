main_Velm_2_4_years <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Velmeshev_2022_2_4_years/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/plot_high_variable_genes_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2_4_years/cell_info_Velmeshev_2022_2_4_years.csv")
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

Velm_2_4_years <- list()
for (v in runs) {
  Velm_2_4_years_input_seurat <- SCENICInput(main_Velm_2_4_years, F, v)
  Velm_2_4_years_seurat_list <- list()
  for (Velm_2_4_years_id in names(Velm_2_4_years_input_seurat)) {
    colnames(Velm_2_4_years_input_seurat[[Velm_2_4_years_id]]) <- str_replace_all(colnames(Velm_2_4_years_input_seurat[[Velm_2_4_years_id]]), "\\.", "-")
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(Velm_2_4_years_input_seurat[[Velm_2_4_years_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    Velm_2_4_years_seurat <- SCENICSeuratPlots(Velm_2_4_years_input_seurat, metadata_id, Velm_2_4_years_id)
    Velm_2_4_years_seurat_list <- append(Velm_2_4_years_seurat_list, list(Velm_2_4_years_seurat))
  }
  names(Velm_2_4_years_seurat_list) <- names(Velm_2_4_years_input_seurat)
  Velm_2_4_years <- append(Velm_2_4_years, list(Velm_2_4_years_seurat_list))
}
names(Velm_2_4_years) <- runs


k_clusters <- list("_1" = c(12, 12),
                   "_2" = c(12, 12),
                   "_3" = c(12, 12))

Velm_2_4_years_final <- list()
for (v in runs) { 
  Velm_2_4_years_seurat_list <- Velm_2_4_years[[v]]
  names(k_clusters[[v]]) <- names(Velm_2_4_years_seurat_list)
  for (Velm_2_4_years_id in names(Velm_2_4_years_seurat_list)) {
    #Velm_2_4_years_final[[v]][[Velm_2_4_years_id]] <- SCENICClustering(main_Velm_2_4_years, F, Velm_2_4_years_seurat_list[[Velm_2_4_years_id]], k_clusters[[v]][[Velm_2_4_years_id]], ct_order)
    # if also need to plot the UMAPs
    Velm_2_4_years_final[[v]][[Velm_2_4_years_id]] <- SCENICClustering(main_Velm_2_4_years, F, Velm_2_4_years_seurat_list[[Velm_2_4_years_id]], k_clusters[[v]][[Velm_2_4_years_id]], ct_order, plot_flag = "yes", Velm_2_4_years_id)
    SCENICMarkers(main_Velm_2_4_years, F, Velm_2_4_years[[v]][[Velm_2_4_years_id]], Velm_2_4_years_id)
  }
}

rm(Velm_2_4_years, Velm_2_4_years_id, Velm_2_4_years_input_seurat, Velm_2_4_years_seurat, Velm_2_4_years_seurat_list, k_clusters, metadata_id)

saveRDS(Velm_2_4_years_final, paste0(main_Velm_2_4_years, "seurat_files.rds"))

Velm_2_4_years_final <- readRDS(paste0(main_Velm_2_4_years, "seurat_files.rds"))

for (run in names(Velm_2_4_years_final)) {
  for (sex in names(Velm_2_4_years_final[[run]])) {
    SCENICUmap(main_Velm_2_4_years, F, Velm_2_4_years_final[[run]][[sex]], ct_order, sex)
  }
}

### Calculates Markers in each SeuratObject

Velm_2_4_years_markers <- SCENICInputMarkers(main_Velm_2_4_years, F, pval_thresh, FC_thresh)
Velm_2_4_years_10 <- SCENICtop10genes(Velm_2_4_years_markers, F)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Velm_2_4_years, F, Velm_2_4_years_final, Velm_2_4_years_10, ct_order)
#HmpSCENICAll(main_Velm_2_4_years, F, Velm_2_4_years_final, Velm_2_4_years_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

Velm_2_4_years_scenic <- ImportSCENICresults(main_Velm_2_4_years, F, "1_GRN", proj_order = "no")
#TfTgSeuratExpression(main_Velm_2_4_years, F, Velm_2_4_years_scenic, Velm_2_4_years_final, ct_order)
TfTgSeuratExpression(main_Velm_2_4_years, F, Velm_2_4_years_scenic, Velm_2_4_years_final, ct_order, 100)

Velm_2_4_years_tf_list <- SCENICExtractGRN(Velm_2_4_years_scenic, F, "TF", 100)
ExtractDiffGRN(main_Velm_2_4_years, F, Velm_2_4_years_tf_list, "TF")
SCENICPlotGRN(main_Velm_2_4_years, F, Velm_2_4_years_tf_list, "TF")

Velm_2_4_years_tg_list <- SCENICExtractGRN(Velm_2_4_years_scenic, F, "target", 50)
ExtractDiffGRN(main_Velm_2_4_years, F, Velm_2_4_years_tg_list, "Target")
SCENICPlotGRN(main_Velm_2_4_years, F, Velm_2_4_years_tg_list, "Target")

#####  TFs and Targets expression in original SeuratObject

Velm_2_4_years <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_2_4_years.rds")
Velm_2_4_years@meta.data$ct_sex <- paste(Velm_2_4_years@meta.data$cluster_final, Velm_2_4_years@meta.data$sex, sep="_")

Velm_2_4_years_tf <- read.csv(paste0(main_Velm_2_4_years, "5_outputs/different_TF_between_sexes.csv"))
RidgeTFTG(main_Velm_2_4_years, Velm_2_4_years, Velm_2_4_years_tf$gene_id, "ct_sex", "TF")
RidgeTFTG(main_Velm_2_4_years, Velm_2_4_years, Velm_2_4_years_tf$gene_id, "sex", "TF")

Velm_2_4_years_tg <- read.csv(paste0(main_Velm_2_4_years, "5_outputs/different_Target_between_sexes.csv"))
RidgeTFTG(main_Velm_2_4_years, Velm_2_4_years, Velm_2_4_years_tg$gene_id, "ct_sex", "Target")
RidgeTFTG(main_Velm_2_4_years, Velm_2_4_years, Velm_2_4_years_tg$gene_id, "sex", "Target")




#####  Regulons

# Velm_2_4_years
Velm_2_4_years_auc <- ImportSCENICresults(main_Velm_2_4_years, F, "3_AUCell", proj_order = "yes")
Velm_2_4_years_reg_list <- SCENICExtractRegulons(Velm_2_4_years_auc, F)
SCENICPlotRegulons(main_Velm_2_4_years, F, Velm_2_4_years_reg_list)

########## Number of TF-TG pairs between F and M of same project
Velm_2_4_years_overlapTFTG <- SCENICOverlapTfTg(Velm_2_4_years_scenic, F, analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_2_4_years, F, Velm_2_4_years_overlapTFTG)
Velm_2_4_years_overlapTFTG <- SCENICOverlapTfTg(Velm_2_4_years_scenic, F, threshold = 10000,  analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_2_4_years, F, Velm_2_4_years_overlapTFTG, threshold = 10000)

#####  Gene Variability

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2_4_years/cell_info_Velmeshev_2022_2_4_years.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)

cell_info[which(cell_info$sex=="Female"), "sex"] <- "F"
cell_info[which(cell_info$sex=="Male"), "sex"] <- "M"

top2000 <- readRDS(paste0(main_Velm_2_4_years, "top_2000_SD_expr_matrix_Velmeshev_2022_2_4_years.rds"))

SexSD(main_Velm_2_4_years, cell_info, top2000)

