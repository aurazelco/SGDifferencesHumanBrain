main_Velm_Adult <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Velmeshev_2022_Adult/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/plot_high_variable_genes_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_Adult/cell_info_Velmeshev_2022_Adult.csv")
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

Velm_Adult <- list()
for (v in runs) {
  Velm_Adult_input_seurat <- SCENICInput(main_Velm_Adult, F, v)
  Velm_Adult_seurat_list <- list()
  for (Velm_Adult_id in names(Velm_Adult_input_seurat)) {
    colnames(Velm_Adult_input_seurat[[Velm_Adult_id]]) <- str_replace_all(colnames(Velm_Adult_input_seurat[[Velm_Adult_id]]), "\\.", "-")
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(Velm_Adult_input_seurat[[Velm_Adult_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    Velm_Adult_seurat <- SCENICSeuratPlots(Velm_Adult_input_seurat, metadata_id, Velm_Adult_id)
    Velm_Adult_seurat_list <- append(Velm_Adult_seurat_list, list(Velm_Adult_seurat))
  }
  names(Velm_Adult_seurat_list) <- names(Velm_Adult_input_seurat)
  Velm_Adult <- append(Velm_Adult, list(Velm_Adult_seurat_list))
}
names(Velm_Adult) <- runs


k_clusters <- list("_1" = c(13, 13),
                   "_2" = c(12, 12),
                   "_3" = c(13, 13))

Velm_Adult_final <- list()
for (v in runs) { 
  Velm_Adult_seurat_list <- Velm_Adult[[v]]
  names(k_clusters[[v]]) <- names(Velm_Adult_seurat_list)
  for (Velm_Adult_id in names(Velm_Adult_seurat_list)) {
    #Velm_Adult_final[[v]][[Velm_Adult_id]] <- SCENICClustering(main_Velm_Adult, F, Velm_Adult_seurat_list[[Velm_Adult_id]], k_clusters[[v]][[Velm_Adult_id]], ct_order)
    # if also need to plot the UMAPs
    Velm_Adult_final[[v]][[Velm_Adult_id]] <- SCENICClustering(main_Velm_Adult, F, Velm_Adult_seurat_list[[Velm_Adult_id]], k_clusters[[v]][[Velm_Adult_id]], ct_order, plot_flag = "yes", Velm_Adult_id)
    SCENICMarkers(main_Velm_Adult, F, Velm_Adult[[v]][[Velm_Adult_id]], Velm_Adult_id)
  }
}

rm(Velm_Adult, Velm_Adult_id, Velm_Adult_input_seurat, Velm_Adult_seurat, Velm_Adult_seurat_list, k_clusters, metadata_id)

saveRDS(Velm_Adult_final, paste0(main_Velm_Adult, "seurat_files.rds"))

Velm_Adult_final <- readRDS(paste0(main_Velm_Adult, "seurat_files.rds"))

for (run in names(Velm_Adult_final)) {
  for (sex in names(Velm_Adult_final[[run]])) {
    SCENICUmap(main_Velm_Adult, F, Velm_Adult_final[[run]][[sex]], ct_order, sex)
  }
}

### Calculates Markers in each SeuratObject

Velm_Adult_markers <- SCENICInputMarkers(main_Velm_Adult, F, pval_thresh, FC_thresh)
Velm_Adult_10 <- SCENICtop10genes(Velm_Adult_markers, F)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_Velm_Adult, F, Velm_Adult_final, Velm_Adult_10, ct_order)
#HmpSCENICAll(main_Velm_Adult, F, Velm_Adult_final, Velm_Adult_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

Velm_Adult_scenic <- ImportSCENICresults(main_Velm_Adult, F, "1_GRN", proj_order = "no")
#TfTgSeuratExpression(main_Velm_Adult, F, Velm_Adult_scenic, Velm_Adult_final, ct_order)
TfTgSeuratExpression(main_Velm_Adult, F, Velm_Adult_scenic, Velm_Adult_final, ct_order, 100)

Velm_Adult_tf_list <- SCENICExtractGRN(Velm_Adult_scenic, F, "TF", 100)
ExtractDiffGRN(main_Velm_Adult, F, Velm_Adult_tf_list, "TF")
SCENICPlotGRN(main_Velm_Adult, F, Velm_Adult_tf_list, "TF")

Velm_Adult_tg_list <- SCENICExtractGRN(Velm_Adult_scenic, F, "target", 50)
ExtractDiffGRN(main_Velm_Adult, F, Velm_Adult_tg_list, "Target")
SCENICPlotGRN(main_Velm_Adult, F, Velm_Adult_tg_list, "Target")

#####  TFs and Targets expression in original SeuratObject

Velm_adult <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_Adult.rds")
Velm_adult@meta.data$ct_sex <- paste(Velm_adult@meta.data$cluster_final, Velm_adult@meta.data$sex, sep="_")

Velm_adult_tf <- read.csv(paste0(main_Velm_Adult, "5_outputs/different_TF_between_sexes.csv"))
RidgeTFTG(main_Velm_Adult, Velm_adult, Velm_adult_tf$gene_id, "ct_sex", "TF")
RidgeTFTG(main_Velm_Adult, Velm_adult, Velm_adult_tf$gene_id, "sex", "TF")


Velm_adult_tg <- read.csv(paste0(main_Velm_Adult, "5_outputs/different_Target_between_sexes.csv"))
RidgeTFTG(main_Velm_Adult, Velm_adult, Velm_adult_tg$gene_id, "ct_sex", "Target")
RidgeTFTG(main_Velm_Adult, Velm_adult, Velm_adult_tg$gene_id, "sex", "Target")


#####  Regulons

# Velm_Adult
Velm_Adult_auc <- ImportSCENICresults(main_Velm_Adult, F, "3_AUCell", proj_order = "yes")
Velm_Adult_reg_list <- SCENICExtractRegulons(Velm_Adult_auc, F)
SCENICPlotRegulons(main_Velm_Adult, F, Velm_Adult_reg_list)

########## Number of TF-TG pairs between F and M of same project
Velm_Adult_overlapTFTG <- SCENICOverlapTfTg(Velm_Adult_scenic, F, analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_Adult, F, Velm_Adult_overlapTFTG)
Velm_Adult_overlapTFTG <- SCENICOverlapTfTg(Velm_Adult_scenic, F, threshold = 10000,  analysis_type = "Velmeshev")
SCENICPlotOverlapTfTg(main_Velm_Adult, F, Velm_Adult_overlapTFTG, threshold = 10000)


#####  Gene Variability

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/SCENIC/plot_high_variable_genes_func.R")

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_Adult/cell_info_Velmeshev_2022_Adult.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("sex", "ct"), sep="_", remove=F)

cell_info[which(cell_info$sex=="Female"), "sex"] <- "F"
cell_info[which(cell_info$sex=="Male"), "sex"] <- "M"

top2000 <- readRDS(paste0(main_Velm_Adult, "top_2000_SD_expr_matrix_Velmeshev_2022_Adult.rds"))

SexSD(main_Velm_Adult, cell_info, top2000)

