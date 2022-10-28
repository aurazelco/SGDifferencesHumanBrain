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
runs <- c( "_1", "_2", "_3")

########## Plot cell types back to umap/tsne

# Normal

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
    # all_norm_final[[v]][[norm_id]] <- SCENICClustering(main, sub_disease[3], norm_seurat_list[[norm_id]], k_clusters[[v]][[norm_id]], ct_order, "yes")
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

#####  TFs and TGs expression in SeuratObjects

norm_scenic <- SCENICresultsSeurat(main, sub_disease[3], "1_GRN")
SCENICTfTg(main, sub_disease[3], norm_scenic, all_norm_final, ct_order)
SCENICTfTg(main, sub_disease[3], norm_scenic, all_norm_final, ct_order, 100)


# AD

#all_ad <- list()
#for (v in runs) {
#  ad_input_seurat <- SCENICInputSeurat(main, sub_disease[1], v)
#  ad_seurat_list <- list()
#  for (ad_id in names(ad_input_seurat)) {
#    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(ad_input_seurat[[ad_id]])),]
#    rownames(metadata_id) <- metadata_id$cell_id
#    metadata_id$cell_id <- NULL
#    ad_seurat <- SCENICSeuratPlots(ad_input_seurat, metadata_id, ad_id)
#    ad_seurat_list <- append(ad_seurat_list, list(ad_seurat))
#  }
#  names(ad_seurat_list) <- names(ad_input_seurat)
#  all_ad <- append(all_ad, list(ad_seurat_list))
#}
#names(all_ad) <- runs

#k_clusters <- list("_1" = c(13, 13, 10, 12),
#                   "_2" = c(13, 13, 12, 11),
#                   "_3" = c(11, 12, 12, 11))

#all_ad_final <- list()
#for (v in runs) { 
#  ad_seurat_list <- all_ad[[v]]
#  names(k_clusters[[v]]) <- names(ad_seurat_list)
#  for (ad_id in names(ad_seurat_list)) {
#    #all_ad_final[[v]][[ad_id]] <- SCENICClustering(main, sub_disease[1], ad_seurat_list[[ad_id]], k_clusters[[v]][[ad_id]], ct_order)
#    # if also need to plot the UMAPs
#    all_ad_final[[v]][[ad_id]] <- SCENICClustering(main, sub_disease[1], ad_seurat_list[[ad_id]], k_clusters[[v]][[ad_id]], ct_order, "yes")
#    SCENICMarkers(main, sub_disease[1], all_ad[[v]][[ad_id]])
#  }
#}
#rm(all_ad)

#saveRDS(all_ad_final, paste0(main, sub_disease[1], "/seurat_files.rds"))
all_ad_final <- readRDS(paste0(main, sub_disease[1], "/seurat_files.rds"))

### Calculates Markers in each SeuratObject

ad_markers <- SCENICInputMarkers(main, sub_disease[1], pval_thresh, FC_thresh)
ad_10 <- SCENICtop10genes(ad_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main, sub_disease[1], all_ad_final, ad_10, ct_order)
HmpSCENICAll(main, sub_disease[1], all_ad_final, ad_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

ad_scenic <- SCENICresultsSeurat(main, sub_disease[1], "1_GRN")
SCENICTfTg(main, sub_disease[1], ad_scenic, all_ad_final, ct_order)
SCENICTfTg(main, sub_disease[1], ad_scenic, all_ad_final, ct_order, 100)


# MS
#all_ms <- list()
#for (v in runs) {
#  ms_input_seurat <- SCENICInputSeurat(main, sub_disease[2], v)
#  ms_seurat_list <- list()
#  for (ms_id in names(ms_input_seurat)) {
#    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(ms_input_seurat[[ms_id]])),]
#    rownames(metadata_id) <- metadata_id$cell_id
#    metadata_id$cell_id <- NULL
#    ms_seurat <- SCENICSeuratPlots(ms_input_seurat, metadata_id, ms_id)
#    ms_seurat_list <- append(ms_seurat_list, list(ms_seurat))
#  }
#  names(ms_seurat_list) <- names(ms_input_seurat)
#  all_ms <- append(all_ms, list(ms_seurat_list))
#}
#names(all_ms) <- runs

#k_clusters <- list("_1" = c(15, 15),
#                   "_2" = c(15, 15),
#                   "_3" = c(15, 15))

#all_ms_final <- list()
#for (v in runs) { 
#  ms_seurat_list <- all_ms[[v]]
#  names(k_clusters[[v]]) <- names(ms_seurat_list)
#  for (ms_id in names(ms_seurat_list)) {
#    #all_ms_final[[v]][[ms_id]] <- SCENICClustering(main, sub_disease[2], ms_seurat_list[[ms_id]], k_clusters[[v]][[ms_id]], ct_order)
#    #if also need to plot the UMAPs
#    all_ms_final[[v]][[ms_id]] <- SCENICClustering(main, sub_disease[2], ms_seurat_list[[ms_id]], k_clusters[[v]][[ms_id]], ct_order, "yes")
#    SCENICMarkers(main, sub_disease[2], all_ms[[v]][[ms_id]])
#  }
#}
#rm(all_ms)

#saveRDS(all_ms_final, paste0(main, sub_disease[2], "/seurat_files.rds"))

all_ms_final <- readRDS(paste0(main, sub_disease[2], "/seurat_files.rds"))

### Calculates Markers in each SeuratObject

ms_markers <- SCENICInputMarkers(main, sub_disease[2], pval_thresh, FC_thresh)
ms_10 <- SCENICtop10genes(ms_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main, sub_disease[2], all_ms_final, ms_10, ct_order)
HmpSCENICAll(main, sub_disease[2], all_ms_final, ms_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

ms_scenic <- SCENICresultsSeurat(main, sub_disease[2], "1_GRN")
SCENICTfTg(main, sub_disease[2], ms_scenic, all_ms_final, ct_order)
SCENICTfTg(main, sub_disease[2], ms_scenic, all_ms_final, ct_order, 100)






#####  Regulons

# check regulons distr
#for (n in 1:nrow(norm_auc[["_1"]][["GSE157827_F_1"]])) {
#  print(hist(as.numeric(norm_auc[["_1"]][["GSE157827_F_1"]][n,2:ncol(norm_auc[["_1"]][["GSE157827_F_1"]])])))
#}

# Normal
norm_auc <- SCENICresultsSeurat(main, sub_disease[3], "3_AUCell")
norm_reg_list <- SCENICExtractRegulons(norm_auc)
SCENICPlotRegulons(main, sub_disease[3], norm_reg_list)

# AD
ad_auc <- SCENICresultsSeurat(main, sub_disease[1], "3_AUCell")
ad_reg_list <- SCENICExtractRegulons(ad_auc)
SCENICPlotRegulons(main, sub_disease[1], ad_reg_list)

# MS
ms_auc <- SCENICresultsSeurat(main, sub_disease[2], "3_AUCell")
ms_reg_list <- SCENICExtractRegulons(ms_auc)
SCENICPlotRegulons(main, sub_disease[2], ms_reg_list)







########## Number of TF-TG pairs between F and M of same project

# Normal
norm_scenic <- SCENICresultsSeurat(main, sub_disease[3], "1_GRN", "yes")
norm_overlapTFTG <- SCENICOverlapTfTg(norm_scenic)
SCENICPlotOverlapTfTg(main, sub_disease[3], norm_overlapTFTG)

# AD
ad_scenic <- SCENICresultsSeurat(main, sub_disease[1], "1_GRN", "yes")
ad_overlapTFTG <- SCENICOverlapTfTg(ad_scenic)
SCENICPlotOverlapTfTg(main, sub_disease[1],ad_overlapTFTG)

# MS
ms_scenic <- SCENICresultsSeurat(main, sub_disease[2], "1_GRN", "yes")
ms_overlapTFTG <- SCENICOverlapTfTg(ms_scenic)
SCENICPlotOverlapTfTg(main, sub_disease[2], ms_overlapTFTG)











###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

########### ONLY TESTED ON NORMAL 1 DF - DO NOT RUN FOR ALL DISEASES

########## Overlap among runs for cts

ct_lists <- c("microglia", "OPC", "oligodendrocyte")
top_lists <- c("no", 1000, 2000, 5000, 10000, 20000)

norm_overlap <- SCENICoverlap(main, sub_disease[3])
PlotOverlapRuns(main, sub_disease[3], norm_overlap, top_lists, ct_lists)






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

