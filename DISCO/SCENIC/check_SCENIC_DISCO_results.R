main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/outputs/"
sub_disease <- list.dirs(main_DISCO, full.names = FALSE, recursive = FALSE)

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/check_SCENIC_results_func.R")
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/SCENIC/plot_high_variable_genes_func.R")

########## Important files

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/extra_files/cell_info.csv")
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

all_norm <- list()
for (v in runs) {
  norm_input_seurat <- SCENICInput(main_DISCO, sub_disease[3], v)
  norm_seurat_list <- list()
  for (norm_id in names(norm_input_seurat)) {
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[norm_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    norm_seurat <- SCENICSeuratPlots(norm_input_seurat, metadata_id, norm_id)
    norm_seurat_list <- append(norm_seurat_list, list(norm_seurat))
  }
  names(norm_seurat_list) <- names(norm_input_seurat)
  all_norm <- append(all_norm, list(norm_seurat_list))
}
names(all_norm) <- runs

k_clusters <- list("_1" = c(15, 15, 10, 10, 10, 12),
                   "_2" = c(12, 10, 10, 10, 12, 15),
                   "_3" = c(15, 10, 10, 10, 12, 12))

all_norm_final <- list()
for (v in runs) { 
  norm_seurat_list <- all_norm[[v]]
  names(k_clusters[[v]]) <- names(norm_seurat_list)
  for (norm_id in names(norm_seurat_list)) {
    all_norm_final[[v]][[norm_id]] <- SCENICClustering(main_DISCO, sub_disease[3], norm_seurat_list[[norm_id]], k_clusters[[v]][[norm_id]], ct_order)
    # if also need to plot the UMAPs
    all_norm_final[[v]][[norm_id]] <- SCENICClustering(main_DISCO, sub_disease[3], norm_seurat_list[[norm_id]], k_clusters[[v]][[norm_id]], ct_order, "yes")
    SCENICMarkers(main_DISCO, sub_disease[3], all_norm[[v]][[norm_id]])
  }
}
rm(all_norm)

saveRDS(all_norm_final, paste0(main_DISCO, sub_disease[3], "/seurat_files.rds"))


### Calculates Markers in each SeuratObject

norm_markers <- SCENICInputMarkers(main_DISCO, sub_disease[3], pval_thresh, FC_thresh)
norm_10 <- SCENICtop10genes(norm_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_DISCO, sub_disease[3], all_norm_final, norm_10, ct_order)
#HmpSCENICAll(main_DISCO, sub_disease[3], all_norm_final, norm_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

norm_scenic <- ImportSCENICresults(main_DISCO, sub_disease[3], "1_GRN")
#TfTgSeuratExpression(main_DISCO, sub_disease[3], norm_scenic, all_norm_final, ct_order)
TfTgSeuratExpression(main_DISCO, sub_disease[3], norm_scenic, all_norm_final, ct_order, 100)

norm_tf_list <- SCENICExtractGRN(norm_scenic, sub_disease[3], "TF", 100)
norm_tf_sex <- ExtractDiffGRN(main_DISCO, sub_disease[3], norm_tf_list, "TF")
SCENICPlotGRN(main_DISCO, sub_disease[3], norm_tf_list, "TF")

norm_tg_list <- SCENICExtractGRN(norm_scenic, sub_disease[3], "target", 50)
norm_tg_sex <- ExtractDiffGRN(main_DISCO, sub_disease[3], norm_tg_list, "Target")
SCENICPlotGRN(main_DISCO, sub_disease[3], norm_tg_list, "Target")


#####  Gene Variability

SexSD(main, cell_info, top2000)


# AD

all_ad <- list()
for (v in runs) {
  ad_input_seurat <- SCENICInput(main_DISCO, sub_disease[1], v)
  ad_seurat_list <- list()
  for (ad_id in names(ad_input_seurat)) {
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(ad_input_seurat[[ad_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    ad_seurat <- SCENICSeuratPlots(ad_input_seurat, metadata_id, ad_id)
    ad_seurat_list <- append(ad_seurat_list, list(ad_seurat))
  }
  names(ad_seurat_list) <- names(ad_input_seurat)
  all_ad <- append(all_ad, list(ad_seurat_list))
}
names(all_ad) <- runs

k_clusters <- list("_1" = c(13, 13, 10, 12),
                   "_2" = c(13, 13, 12, 11),
                   "_3" = c(11, 12, 12, 11))

all_ad_final <- list()
for (v in runs) { 
  ad_seurat_list <- all_ad[[v]]
  names(k_clusters[[v]]) <- names(ad_seurat_list)
  for (ad_id in names(ad_seurat_list)) {
    all_ad_final[[v]][[ad_id]] <- SCENICClustering(main_DISCO, sub_disease[1], ad_seurat_list[[ad_id]], k_clusters[[v]][[ad_id]], ct_order)
    # if also need to plot the UMAPs
    all_ad_final[[v]][[ad_id]] <- SCENICClustering(main_DISCO, sub_disease[1], ad_seurat_list[[ad_id]], k_clusters[[v]][[ad_id]], ct_order, "yes")
    SCENICMarkers(main_DISCO, sub_disease[1], all_ad[[v]][[ad_id]])
  }
}
rm(all_ad)

saveRDS(all_ad_final, paste0(main_DISCO, sub_disease[1], "/seurat_files.rds"))
all_ad_final <- readRDS(paste0(main_DISCO, sub_disease[1], "/seurat_files.rds"))

### Calculates Markers in each SeuratObject

ad_markers <- SCENICInputMarkers(main_DISCO, sub_disease[1], pval_thresh, FC_thresh)
ad_10 <- SCENICtop10genes(ad_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_DISCO, sub_disease[1], all_ad_final, ad_10, ct_order)
#HmpSCENICAll(main_DISCO, sub_disease[1], all_ad_final, ad_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

ad_scenic <- ImportSCENICresults(main_DISCO, sub_disease[1], "1_GRN")
#TfTgSeuratExpression(main_DISCO, sub_disease[1], ad_scenic, all_ad_final, ct_order)
TfTgSeuratExpression(main_DISCO, sub_disease[1], ad_scenic, all_ad_final, ct_order, 100)

ad_tf_list <- SCENICExtractGRN(ad_scenic, sub_disease[1], "TF", 100)
ExtractDiffGRN(main_DISCO, sub_disease[1], ad_tf_list, "TF")
SCENICPlotGRN(main_DISCO, sub_disease[1], ad_tf_list, "TF")

ad_tg_list <- SCENICExtractGRN(ad_scenic, sub_disease[1], "target", 50)
ExtractDiffGRN(main_DISCO, sub_disease[1], ad_tg_list, "Target")
SCENICPlotGRN(main_DISCO, sub_disease[1], ad_tg_list, "Target")


# MS
all_ms <- list()
for (v in runs) {
  ms_input_seurat <- SCENICInput(main_DISCO, sub_disease[2], v)
  ms_seurat_list <- list()
  for (ms_id in names(ms_input_seurat)) {
    metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(ms_input_seurat[[ms_id]])),]
    rownames(metadata_id) <- metadata_id$cell_id
    metadata_id$cell_id <- NULL
    ms_seurat <- SCENICSeuratPlots(ms_input_seurat, metadata_id, ms_id)
    ms_seurat_list <- append(ms_seurat_list, list(ms_seurat))
  }
  names(ms_seurat_list) <- names(ms_input_seurat)
  all_ms <- append(all_ms, list(ms_seurat_list))
}
names(all_ms) <- runs

k_clusters <- list("_1" = c(15, 15),
                   "_2" = c(15, 15),
                   "_3" = c(15, 15))

all_ms_final <- list()
for (v in runs) { 
  ms_seurat_list <- all_ms[[v]]
  names(k_clusters[[v]]) <- names(ms_seurat_list)
  for (ms_id in names(ms_seurat_list)) {
    all_ms_final[[v]][[ms_id]] <- SCENICClustering(main_DISCO, sub_disease[2], ms_seurat_list[[ms_id]], k_clusters[[v]][[ms_id]], ct_order)
    #if also need to plot the UMAPs
    all_ms_final[[v]][[ms_id]] <- SCENICClustering(main_DISCO, sub_disease[2], ms_seurat_list[[ms_id]], k_clusters[[v]][[ms_id]], ct_order, "yes")
    SCENICMarkers(main_DISCO, sub_disease[2], all_ms[[v]][[ms_id]])
  }
}
rm(all_ms)

saveRDS(all_ms_final, paste0(main_DISCO, sub_disease[2], "/seurat_files.rds"))

all_ms_final <- readRDS(paste0(main_DISCO, sub_disease[2], "/seurat_files.rds"))

### Calculates Markers in each SeuratObject
ms_markers <- SCENICInputMarkers(main_DISCO, sub_disease[2], pval_thresh, FC_thresh)
ms_10 <- SCENICtop10genes(ms_markers)

### Calculates Markers in each SeuratObject

HmpSCENIC(main_DISCO, sub_disease[2], all_ms_final, ms_10, ct_order)
#HmpSCENICAll(main_DISCO, sub_disease[2], all_ms_final, ms_markers, ct_order)

#####  TFs and TGs expression in SeuratObjects

ms_scenic <- ImportSCENICresults(main_DISCO, sub_disease[2], "1_GRN")
#TfTgSeuratExpression(main_DISCO, sub_disease[2], ms_scenic, all_ms_final, ct_order)
TfTgSeuratExpression(main_DISCO, sub_disease[2], ms_scenic, all_ms_final, ct_order, 100)


ms_tf_list <- SCENICExtractGRN(ms_scenic, sub_disease[2], "TF", 100)
ExtractDiffGRN(main_DISCO, sub_disease[2], ms_tf_list, "TF")
SCENICPlotGRN(main_DISCO, sub_disease[2], ms_tf_list, "TF")

ms_tg_list <- SCENICExtractGRN(ms_scenic, sub_disease[2], "target", 50)
ExtractDiffGRN(main_DISCO, sub_disease[2], ms_tg_list, "Target")
SCENICPlotGRN(main_DISCO, sub_disease[2], ms_tg_list, "Target")




#####  Regulons

# Normal
norm_auc <- ImportSCENICresults(main_DISCO, sub_disease[3], "3_AUCell")
norm_reg_list <- SCENICExtractRegulons(norm_auc)
SCENICPlotRegulons(main_DISCO, sub_disease[3], norm_reg_list)

# AD
ad_auc <- ImportSCENICresults(main_DISCO, sub_disease[1], "3_AUCell")
ad_reg_list <- SCENICExtractRegulons(ad_auc)
SCENICPlotRegulons(main_DISCO, sub_disease[1], ad_reg_list)

# MS
ms_auc <- ImportSCENICresults(main_DISCO, sub_disease[2], "3_AUCell")
ms_reg_list <- SCENICExtractRegulons(ms_auc)
SCENICPlotRegulons(main_DISCO, sub_disease[2], ms_reg_list)



########## Number of TF-TG pairs between F and M of same project

# Normal
norm_scenic <- ImportSCENICresults(main_DISCO, sub_disease[3], "1_GRN", "yes")
norm_overlapTFTG <- SCENICOverlapTfTg(norm_scenic, sub_disease[3])
SCENICPlotOverlapTfTg(main_DISCO, sub_disease[3], norm_overlapTFTG)
norm_overlapTFTG <- SCENICOverlapTfTg(norm_scenic, sub_disease[3], threshold = 10000)
SCENICPlotOverlapTfTg(main_DISCO, sub_disease[3], norm_overlapTFTG, threshold = 10000)

# AD
ad_scenic <- ImportSCENICresults(main_DISCO, sub_disease[1], "1_GRN", "yes")
ad_overlapTFTG <- SCENICOverlapTfTg(ad_scenic, sub_disease[1], threshold = 10000)
SCENICPlotOverlapTfTg(main_DISCO, sub_disease[1], ad_overlapTFTG, threshold = 10000)

# MS
ms_scenic <- ImportSCENICresults(main_DISCO, sub_disease[2], "1_GRN", "yes")
ms_overlapTFTG <- SCENICOverlapTfTg(ms_scenic, sub_disease[2], threshold = 10000)
SCENICPlotOverlapTfTg(main_DISCO, sub_disease[2], ms_overlapTFTG, threshold = 10000)


#####  Gene Variability - ON KJEMPEFURU

#####  Gene Variability

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/SCENIC/plot_high_variable_genes_func.R")

main_scenic <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/outputs/"

cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs/cell_info.csv")
cell_info$X <- NULL
cell_info <- separate(cell_info, og_group, into=c("proj", "sex","disease", "ct"), sep="_", remove=F)
library(stringr)
cell_info$ct <- str_replace_all(cell_info$ct, c("/"="_", " "="_"))

sub_disease <- list.dirs(main_scenic, full.names = F, recursive = F)

for (dis_type in sub_disease) {
  print(dis_type)
  main_scenic_dis <- paste0(main_scenic, dis_type, "/")
  dir.create(main_scenic_dis, recursive = T, showWarnings = F)
  top2000_dis <- readRDS(paste0(main_scenic_dis, "top_2000_SD_expr_matrix.rds"))
  cell_info_dis <- subset(cell_info, disease==dis_type)
  if (length(unique(cell_info_dis$proj)) == 1) {
    print("only one project")
    SexSD(main_scenic_dis, cell_info_dis, top2000_dis)
  } else {
    for (proj_id in unique(cell_info_dis$proj)) {
      print(paste0("multiple projects: ", proj_id))
      cell_info_ds_proj <- subset(cell_info_dis, proj==proj_id)
      main_scenic_dis_proj <- paste0(main_scenic_dis, "plots/", proj_id, "/")
      dir.create(main_scenic_dis_proj, recursive = T, showWarnings = F)
      SexSD(main_scenic_dis_proj, cell_info_ds_proj, top2000_dis)
    }
  }
  rm(top2000_dis, cell_info_dis)
}


#####  TFs and Targets expression in original SeuratObject

disco_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/"
disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))


disco_filt.list <- SplitObject(disco_filt, split.by = "disease")
rm(disco_filt)

norm <- disco_filt.list[[1]]
ad <- disco_filt.list[[2]]
ms <- disco_filt.list[[3]]
rm(disco_filt.list)

###### NORMAL

norm@meta.data$ct_sex <- paste(norm@meta.data$ct, norm@meta.data$gender, sep="_")

norm_tf <- read.csv(paste0(main_DISCO, sub_disease[3], "/5_outputs/different_TF_between_sexes.csv"))
norm_tg <- read.csv(paste0(main_DISCO, sub_disease[3], "/5_outputs/different_Target_between_sexes.csv"))

if (length(unique(norm_tf$projs))>1) {
  norm_projs <- SplitObject(norm, split.by = "project_id")
  for (proj in names(norm_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm_projs[[proj]], norm_tf[which(norm_tf$projs==proj), "gene_id"], "ct_sex", paste("TF", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm_projs[[proj]], norm_tf[which(norm_tf$projs==proj), "gene_id"], "gender", paste("TF", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm, norm_tf$gene_id, "ct_sex", "TF")
  RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm, norm_tf$gene_id, "gender", "TF")
}

if (length(unique(norm_tg$projs))>1) {
  norm_projs <- SplitObject(norm, split.by = "project_id")
  for (proj in names(norm_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm_projs[[proj]], norm_tg[which(norm_tg$projs==proj), "gene_id"], "ct_sex", paste("Target", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm_projs[[proj]], norm_tg[which(norm_tg$projs==proj), "gene_id"], "gender", paste("Target", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm, norm_tg$gene_id, "ct_sex", "Target")
  RidgeTFTG(paste0(main_DISCO, sub_disease[3], "/"), norm, norm_tg$gene_id, "gender", "Target")
}


###### AD

ad@meta.data$ct_sex <- paste(ad@meta.data$ct, ad@meta.data$gender, sep="_")

ad_tf <- read.csv(paste0(main_DISCO, sub_disease[1], "/5_outputs/different_TF_between_sexes.csv"))
ad_tg <- read.csv(paste0(main_DISCO, sub_disease[1], "/5_outputs/different_Target_between_sexes.csv"))

if (length(unique(ad_tf$projs))>1) {
  ad_projs <- SplitObject(ad, split.by = "project_id")
  for (proj in names(ad_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad_projs[[proj]], ad_tf[which(ad_tf$projs==proj), "gene_id"], "ct_sex", paste("TF", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad_projs[[proj]], ad_tf[which(ad_tf$projs==proj), "gene_id"], "gender", paste("TF", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad, ad_tf$gene_id, "ct_sex", "TF")
  RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad, ad_tf$gene_id, "gender", "TF")
}

if (length(unique(ad_tg$projs))>1) {
  ad_projs <- SplitObject(ad, split.by = "project_id")
  for (proj in names(ad_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad_projs[[proj]], ad_tg[which(ad_tg$projs==proj), "gene_id"], "ct_sex", paste("Target", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad_projs[[proj]], ad_tg[which(ad_tg$projs==proj), "gene_id"], "gender", paste("Target", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad, ad_tg$gene_id, "ct_sex", "Target")
  RidgeTFTG(paste0(main_DISCO, sub_disease[1], "/"), ad, ad_tg$gene_id, "gender", "Target")
}



###### MS

ms@meta.data$ct_sex <- paste(ms@meta.data$ct, ms@meta.data$gender, sep="_")

ms_tf <- read.csv(paste0(main_DISCO, sub_disease[2], "/5_outputs/different_TF_between_sexes.csv"))
ms_tg <- read.csv(paste0(main_DISCO, sub_disease[2], "/5_outputs/different_Target_between_sexes.csv"))

if (length(unique(ms_tf$projs))>1) {
  ms_projs <- SplitObject(ms, split.by = "project_id")
  for (proj in names(ms_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms_projs[[proj]], ms_tf[which(ms_tf$projs==proj), "gene_id"], "ct_sex", paste("TF", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms_projs[[proj]], ms_tf[which(ms_tf$projs==proj), "gene_id"], "gender", paste("TF", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms, ms_tf$gene_id, "ct_sex", "TF")
  RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms, ms_tf$gene_id, "gender", "TF")
}

if (length(unique(ms_tg$projs))>1) {
  ms_projs <- SplitObject(ms, split.by = "project_id")
  for (proj in names(ms_projs)) {
    print(proj)
    RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms_projs[[proj]], ms_tg[which(ms_tg$projs==proj), "gene_id"], "ct_sex", paste("Target", proj, sep = "_"))
    RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms_projs[[proj]], ms_tg[which(ms_tg$projs==proj), "gene_id"], "gender", paste("Target", proj, sep = "_"))
  }
} else {
  RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms, ms_tg$gene_id, "ct_sex", "Target")
  RidgeTFTG(paste0(main_DISCO, sub_disease[2], "/"), ms, ms_tg$gene_id, "gender", "Target")
}








###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################

########### ONLY TESTED ON NORMAL 1 DF - DO NOT RUN FOR ALL DISEASES

########## Overlap among runs for cts

ct_lists <- c("microglia", "OPC", "oligodendrocyte")
top_lists <- c("no", 1000, 2000, 5000, 10000, 20000)

norm_overlap <- SCENICoverlap(main_DISCO, sub_disease[3])
PlotOverlapRuns(main_DISCO, sub_disease[3], norm_overlap, top_lists, ct_lists)

########## Overlap with scGRNom

scGRNom_sheets <- c("Mic_GRN_with_openchrom", 
                    "Oli_GRN_with_openchrom", 
                    "Mic_GRN_without_openchrom", 
                    "Oli_GRN_without_openchrom")
scGRNom <- list()
for (sheet in scGRNom_sheets) {
  scGRNom_sub <- as.data.frame(read_xlsx(("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/extra_files/scGRNom_suppl_file_2.xlsx"),
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

PlotscGRNomOverlap(main_DISCO, sub_disease[3], norm_overlap, scGRNom_top, scGRNom_ct, flag_scGRNom)

