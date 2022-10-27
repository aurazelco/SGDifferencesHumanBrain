main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/outputs/"
sub_disease <- list.dirs(main, full.names = FALSE, recursive = FALSE)

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/extra_files/cell_info.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("proj", "sex", "disease", "ct"), sep="_", remove=F)
cell_info$ct <- str_replace_all(cell_info$ct, "/", "_")

ct_order <- c(
  "L2_3 EN",               
  "L4 EN",   
  "PLCH1 L4_5 EN", 
  "TSHZ2 L4_5 EN", 
  "L5 EN",       
  "L5_6 EN",       
  "L5b EN",     
  "L6 EN",     
  "pyramidal neuron", 
  "CXCL14 IN",  
  "PVALB IN",                    
  "SST IN",
  "SV2C IN",               
  "VIP IN",  
  "EC", 
  "fibrous astrocyte",
  "protoplasmic astrocyte",
  "OPC", 
  "oligodendrocyte",           
  "microglia"
)

pval_thresh <- 0.05
FC_thresh <- 1.2

########## Plot cell types back to umap/tsne

##### Solution 1 - create a Seurat object for each expression matrix

#runs <- c( "_1", "_2", "_3")

#all_norm <- list()
#for (v in runs) {
#  norm_input_seurat <- SCENICInputSeurat(main, sub_disease[3], v)
#  norm_seurat_list <- list()
#  for (norm_id in names(norm_input_seurat)) {
#    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[norm_id]])),]
#    rownames(metadata_id) <- metadata_id$cell_id
#    metadata_id$cell_id <- NULL
#    norm_seurat <- SCENICSeuratPlots(norm_input_seurat, metadata_id, norm_id)
#    norm_seurat_list <- append(norm_seurat_list, list(norm_seurat))
#  }
#  names(norm_seurat_list) <- names(norm_input_seurat)
#  all_norm <- append(all_norm, list(norm_seurat_list))
#}
#names(all_norm) <- runs

#k_clusters <- list("_1" = c(15, 15, 10, 10, 10, 12),
#                   "_2" = c(12, 10, 10, 10, 12, 15),
#                   "_3" = c(15, 10, 10, 10, 12, 12))

#all_norm_final <- list()
#for (v in runs) { 
#  norm_seurat_list <- all_norm[[v]]
#  names(k_clusters[[v]]) <- names(norm_seurat_list)
#  for (norm_id in names(norm_seurat_list)) {
#    all_norm_final[[v]][[norm_id]] <- SCENICClustering(main, sub_disease[3], norm_seurat_list[[norm_id]], k_clusters[[v]][[norm_id]], ct_order)
    # if also need to plot the UMAPs
    # all_norm2[[v]][[norm_id]] <- SCENICClustering(main, sub_disease[3], norm_seurat_list[[norm_id]], k_clusters[[v]][[norm_id]], ct_order, "yes")
  #  SCENICMarkers(main, sub_disease[3], all_norm[[v]][[norm_id]])
#  }
#}
#rm(all_norm)

#saveRDS(all_norm_final, paste0(main, sub_disease[3], "/seurat_files.rds"))

all_norm_final <- readRDS(paste0(main, sub_disease[3], "/seurat_files.rds"))

### Calculates Markers in each SeuratObject

norm_markers <- SCENICInputMarkers(main, sub_disease[3], pval_thresh, FC_thresh)

norm_10 <- SCENICtop10genes(norm_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main, sub_disease[3], all_norm_final, norm_10, ct_order)

HmpSCENICAll(main, sub_disease[3], all_norm_final, norm_markers, ct_order)


##### Solution 2 - subset from original disco 

#meta_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[1]])),]
#rownames(meta_id) <- meta_id$cell_id
#meta_id$cell_id <- NULL


disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

disco_filt$cell_id <- rownames(disco_filt@meta.data)
disco_filt@meta.data$ct <- str_replace_all(disco_filt@meta.data$ct, "/", "_")
  
for (norm_id in names(norm_input_seurat)) {
  metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[norm_id]])),]
  rownames(metadata_id) <- metadata_id$cell_id
  metadata_id$cell_id <- NULL
  id_sub <- subset(disco_filt, cell_id %in% rownames(metadata_id))
  umap_id_sub <- DimPlot(id_sub, reduction = "umap", group.by = "ct")
  umap_order <- ct_order[which(ct_order %in% levels(umap_id_sub$data$ct))]
  umap_id_sub$data$ct <- factor(umap_id_sub$data$ct, umap_order)
  umap_id_sub$data <- umap_id_sub$data[order(umap_id_sub$data$ct), ]
  pdf(paste0(main, sub_disease[3], "/3_plots/UMAP_", norm_id, "_disco_subset.pdf"))
  print(umap_id_sub  + labs(title = norm_id))
  dev.off()
}

#####  TFs and TGs expression in SeuratObjects

norm_scenic <- SCENICresultsSeurat(main, sub_disease[3], "1_GRN")

SCENICTfTg(main, sub_disease[3], norm_scenic, all_norm_final, ct_order)
SCENICTfTg(main, sub_disease[3], norm_scenic, all_norm_final, ct_order, 100)

########## Overlap among runs for cts
ct_lists <- c("microglia", "OPC", "oligodendrocyte")
top_lists <- c("no", 1000, 2000, 5000, 10000, 20000)

norm_overlap <- SCENICoverlap(main, sub_disease[3])
PlotOverlapRuns(main, sub_disease[3], norm_overlap, top_lists, ct_lists)


#####  Regulons in SeuratObjects

norm_auc <- SCENICresultsSeurat(main, sub_disease[3], "3_AUCell")

for (n in 1:nrow(norm_auc[["_1"]][["GSE157827_F_1"]])) {
  print(hist(as.numeric(norm_auc[["_1"]][["GSE157827_F_1"]][n,2:ncol(norm_auc[["_1"]][["GSE157827_F_1"]])])))
}



########## Overlap with scGRNom

scGRNom_sheets <- c("Mic_GRN_with_openchrom", 
                    "Oli_GRN_with_openchrom", 
                    "Mic_GRN_without_openchrom", 
                    "Oli_GRN_without_openchrom")
scGRNom <- list()
for (sheet in scGRNom_sheets) {
  scGRNom_sub <- as.data.frame(read_xlsx(("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/extra_files/scGRNom_suppl_file_2.xlsx"),
                           sheet = sheet, skip = 1))
  scGRNom <- append(scGRNom, list(scGRNom_sub))
}
names(scGRNom) <- scGRNom_sheets
for (df in names(scGRNom)) {
  scGRNom[[df]]$TF_TG <- paste(scGRNom[[df]]$TF, scGRNom[[df]]$TG, sep="_")
  scGRNom[[df]]$TF_TG <- as.factor(scGRNom[[df]]$TF_TG)
}
scGRNom_ct <- list("microglia"="Mic", "oligodendrocyte"="Oli")
scGRNom_top <- lapply(1:length(names(scGRNom)), function(x) nrow(scGRNom[[x]]))
names(scGRNom_top) <- names(scGRNom)

flag_scGRNom <- c("no", "yes")

PlotscGRNomOverlap(main, sub_disease[3], norm_overlap, scGRNom_top, scGRNom_ct, flag_scGRNom)

########## Heatmap expression of SCENIC TFs and TGs 

#expr_mat_all <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/top_2000_SD_expr_matrix.rds")
#cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/extra_files/cell_info.csv")
#cell_info$X <- NULL
#cell_info <- separate(cell_info, og_group, into=c("proj", "sex", "disease", "ct"), sep="_", remove=F)
#cell_info$ct <- str_replace_all(cell_info$ct, "/", "_")

top_TFs <- c(100, 200)

norm_sort <- SCENICresults(main, sub_disease[3])
norm_input <- SCENICInput(main, sub_disease[3], cell_info)
PlotTfTg(main, sub_disease[3], norm_sort, norm_input, ct_order, top_TFs)
