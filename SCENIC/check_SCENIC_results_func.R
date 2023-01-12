# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Documentation abbreviations:
  # F and M: females and males
  # ct: celltype
  # df: dataframe
  # hmps: heatmaps
  # TF: transcription factor
  # Tg/TG: target

# OBS: this script is sourced in the several check_SCENIC_*_results.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(ggplot2) # to plot
library(readxl) # to read xlsx files
library(stringr) # to modfiy strings
library(stringi) # to modfiy strings
library(tidyr) # to tidy up dfs
library(reshape2) # to re-organize dfs
library(Seurat) # to create a single-cell object
library(SeuratObject) # to create a single-cell object
library(ggpubr) # combine multiple plots into a figure
library(S4Vectors) # 
library(purrr) # 
library(dplyr) # to modify tables
library(matrixStats) # to calculate rapidly some metrics in a matrix

################################## Check randomly sampled SCENIC inputs

# 1. Imports expression matrices used as input for SCENIC, to check if the cells were distanced enough to be separating even with only 100 in each celltype
  # Input: main directory where SCENIC inputs and outputs are found, if there is any sub-folder in the file tree, and which of the three runs we are importing
  # Return: a list of F and M expression matrices from a single sampling to be used to build the SeuratObjects

SCENICInput <- function(main_dir, dis_type, run_v) {
  if (dis_type!=F) {
    input_dfs_path <-  paste0(main_dir, dis_type, "/0_input_dfs/sampled_100_cells_all")
  } else {
    input_dfs_path <-  paste0(main_dir, "0_input_dfs/sampled_100_cells_all")
  }
  input_dfs_files <- list.files(path = input_dfs_path, full.names = T)
  input_dfs <- lapply(input_dfs_files,function(x) {
    read.csv(file = x, 
               sep = ',', 
               header = TRUE)
  })
  names(input_dfs) <- list.files(path = input_dfs_path, full.names = F)
  names(input_dfs) <- str_remove_all(names(input_dfs), ".csv")
  only_1 <- names(input_dfs)[which(grepl(run_v, names(input_dfs)))]
  input_dfs <- input_dfs[only_1]
  for (id in names(input_dfs)) {
    rownames(input_dfs[[id]]) <- input_dfs[[id]]$Genes
    input_dfs[[id]]$Genes <- NULL
  }
  return(input_dfs)
}

# 2. Runs the Seurat clustering from expression matrix (follows tutorial)
  # Input: the list fo expression matrices created before, the corresponding metadata and which one fo the elements in the list to be created 
  # Return: the SeuratObject, ready ot be clustered once the plots generated in this function are evaluated

SCENICSeuratPlots <- function(expr_matrix_list, metadata_id, id) {
  id_seurat <- CreateSeuratObject(as.matrix(expr_matrix_list[[id]]), project = id,
                                    assay = "RNA", meta.data = metadata_id)
  id_seurat[["percent.mt"]] <- PercentageFeatureSet(id_seurat, pattern = "^MT-")
  id_seurat <- subset(id_seurat, subset = percent.mt < 5)
  id_seurat <- NormalizeData(id_seurat)
  id_seurat <- FindVariableFeatures(id_seurat, selection.method = "vst", nfeatures = 2000)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(id_seurat), 10)
  all.genes <- rownames(id_seurat)
  id_seurat <- ScaleData(id_seurat, features = all.genes)
  id_seurat <- RunPCA(id_seurat, features = VariableFeatures(object = id_seurat))
  # Examine and visualize PCA results a few different ways
  print(id_seurat[["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(id_seurat, dims = 1:2, reduction = "pca")
  p1 <- DimHeatmap(id_seurat, dims = 1:15, cells = 500, balanced = TRUE)
  id_seurat <- JackStraw(id_seurat, num.replicate = 100)
  id_seurat <- ScoreJackStraw(id_seurat, dims = 1:20)
  p2 <- JackStrawPlot(id_seurat, dims = 1:15)
  p3 <- ElbowPlot(id_seurat)
  plot<- ggarrange(p2, p3, labels = c("JackStrawPlot", "ElbowPlot"), ncol = 2)
  print(
    annotate_figure(plot, top = text_grob(id, 
                                        color = "black", face = "bold"))
  )
  return(id_seurat)
}

# 3. Runs UMAP on the SeuratObject to define cell clusters
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the SeuratObject, the number of dimensions to use to run th UMAP clustering,
    # the order to display the celltypes, if we have to save the UMAP (see following function), and the UMAP plot filename
  # Return: the SeuratObject clustered and ready for analysis

SCENICClustering <- function(main_dir, dis_type, id_seurat, cluster_num, ct_ordered, plot_flag = "no", file_out_name) {
  id_seurat <- FindNeighbors(id_seurat, dims = 1:cluster_num)
  id_seurat <- FindClusters(id_seurat, resolution = 0.5)
  id_seurat <- RunUMAP(id_seurat, dims = 1:cluster_num)
  if (plot_flag == "yes") {
    SCENICUmap(main_dir, dis_type, id_seurat, ct_ordered, file_out_name)
  }
  return(id_seurat)
}

# 4. Plots the UMAP with the celltype annotation
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the SeuratObject, 
    # the order to display the celltypes, if we have to save the UMAP (see following function), and the UMAP plot filename
  # Return: nothing, saves the plot instead

SCENICUmap <- function (main_dir, dis_type, id_seurat, ct_ordered, file_out_name) {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/Seurat/")
  } else {
    plot_path <- paste0(main_dir, "plots/Seurat/")
  }
  umap_id <- DimPlot(id_seurat, reduction = "umap", group.by = "ct")
  umap_order <- ct_ordered[which(ct_ordered %in% levels(umap_id$data$ct))]
  umap_id$data$ct <- factor(umap_id$data$ct, umap_order)
  umap_id$data <- umap_id$data[order(umap_id$data$ct), ]
  dir.create(plot_path, showWarnings = F, recursive = T)
  pdf(paste0(plot_path,"UMAP_", file_out_name, ".pdf"))
  print(umap_id  + labs(title = file_out_name))
  dev.off()
}

# 5. Finds the markers genes for the celltypes, to double-check that the celltypes are indeed distinct even if 100 cells/ct are used to build the SeuratObject
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the SeuratObject, and the output filename
  # Return: nothing, saves the CSV files instead

SCENICMarkers <- function (main_dir, dis_type, id_seurat, file_out_name) {
  Idents(id_seurat) <- "ct"
  ct_markers <- FindAllMarkers(id_seurat, 
                                    logfc.threshold = 0.25,
                                    min.pct = 0.1)
  if (dis_type!=F) {
    csv_path <- paste0(main_dir, dis_type, "/4_Markers/")
  } else {
    csv_path <- paste0(main_dir, "4_Markers/")
  }
  dir.create(csv_path, showWarnings = F, recursive = T)
  write.csv(ct_markers, file = paste0(csv_path, file_out_name, "_markers.csv"),
            row.names = TRUE)
}

# 6. Filters the genes based on  p-value and fold-change thresholds
  # Input: gene marker df, the p-value and FC thresholds
  # Return: the filtered df

Filter_gene <- function(order.gene.df, pval, FC) {
  logFC <- log2(1/FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  return(gene.sig)
}

# 7. Imports the CSV gene marker files, in order to be plotted into an hmp
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the threshold for the p-value and the fold-change
  # Return: list of gene marker CSV files

SCENICInputMarkers <- function(main_dir, dis_type, pval, FC) {
  if (dis_type!=F) {
    input_markers_path <-  paste0(main_dir, dis_type, "/4_Markers/")
  } else {
    input_markers_path <-  paste0(main_dir, "4_Markers/")
  }
  input_markers_files <- list.files(path = input_markers_path, full.names = T)
  input_markers <- lapply(input_markers_files,function(x) {
    read.table(file = x, 
               sep = ',', 
               header = TRUE)
  })
  input_markers <- lapply(1:length(input_markers), function(x) input_markers[[x]] <- subset(input_markers[[x]], select=-c(1)))
  input_markers <- lapply(1:length(input_markers), function(x) input_markers[[x]] <- Filter_gene(input_markers[[x]], pval, FC))
  names(input_markers) <- list.files(path = input_markers_path, full.names = F)
  names(input_markers) <- str_remove_all(names(input_markers), "_markers.csv")
  for (df in names(input_markers)) {
    input_markers[[df]]$cluster <- as.factor(input_markers[[df]]$cluster)
  }
  runs <- c("_1", "_2", "_3")
  markers_v <- list()
  for (v in runs) {
    v_dfs <- names(input_markers)[which(grepl(v, names(input_markers)))]
    v_list <- input_markers[v_dfs]
    names(v_list) <- v_dfs
    markers_v <- append(markers_v, list(v_list))
  }
  names(markers_v) <- runs
  return(markers_v)
}

# 8. Generate a new table with only the top 10 marker genes for each ct
  # Input: the markers for a specific SeuratObject, if there is any sub-folder in the file tree (if so, there is additional information in the metadata)
  # Return: a df with the top 10 genes for each ct

SCENICtop10genes <- function(input_markers, dis_type) {
  top10 <- data.frame()
  for (v in names(input_markers)) {
    for (df in names(input_markers[[v]])) {
      for (ct in levels(input_markers[[v]][[df]]$cluster)) {
        sub_ct <- subset(input_markers[[v]][[df]], cluster==ct)
        sub_ct <- sub_ct[order(-sub_ct$avg_log2FC),]
        top10 <- rbind(top10, data.frame(rep(paste(df, ct, sep="_"), 10), sub_ct[1:10, "gene"]))
      }
    }
  }
  colnames(top10) <- c("og", "genes")
  for (i in 1:3) {top10$og <- str_replace(top10$og, "_", "/")}
  if (dis_type!=F) {
    top10 <- separate(top10, og, into = c("proj", "sex", "v", "ct"), sep="/", remove=F)
    top10$files <- paste(top10$proj, top10$sex, top10$v, sep = "_")
    col_names <- c("proj", "sex", "v", "ct")
  } else {
    top10 <- separate(top10, og, into = c("sex", "v", "ct"), sep="/", remove=F)
    top10$files <- paste(top10$sex, top10$v, sep = "_")
    col_names <- c("sex", "v", "ct")
  }
  top10[col_names] <- lapply(top10[col_names], as.factor)
  return(top10)
}


# 9. Generates a heatmap displaying the top 10 markers for each ct
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the SeuratObject, the top10 df, the order of the ct, 
    # and if to follow the order or not (sometimes the SeuratObject has issues with the ct order)
  # Return: nothing, saves the plot instead

HmpSCENIC <- function(main_dir, dis_type, input_seurat, top10, ct_ordered, ord_levels="yes") {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/Seurat/")
  } else {
    plot_path <- paste0(main_dir, "plots/Seurat/")
  }
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (run_v in names(input_seurat)) {
    for (file in names(input_seurat[[run_v]])) {
      topgenes <-(top10[which(top10$files==file),])
      topgenes_order <- ct_ordered[which(ct_ordered %in% unique(topgenes$ct))]
      topgenes$ct <- factor(topgenes$ct,topgenes_order)
      topgenes <- topgenes[order(topgenes$ct), ]
      hmp_top_order <- ct_ordered[which(ct_ordered %in% unique(topgenes$ct))]
      Idents(input_seurat[[run_v]][[file]]) <- "ct"
      if (ord_levels=="yes") {
        levels(input_seurat[[run_v]][[file]]) <- hmp_top_order
      }
      hmp_top <- DoHeatmap(input_seurat[[run_v]][[file]], features = topgenes$genes, group.by = "ident", angle = 90, size = 3)
      pdf(paste0(plot_path, "hmp_top_10_", file, ".pdf"), height = 15)
      print(hmp_top  + labs(title = file))
      dev.off()
    }
  }
}

# 10. Generates a heatmap displaying all the markers for each ct
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, the SeuratObject, the markers df, the order of the ct, 
    # and if to follow the order or not (sometimes the SeuratObject has issues with the ct order)
  # Return: nothing, saves the plot instead

HmpSCENICAll <- function(main_dir, dis_type, input_seurat, markers, ct_ordered, ord_levels="yes") {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/Seurat/")
  } else {
    plot_path <- paste0(main_dir, "plots/Seurat/")
  }
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (run_v in names(input_seurat)) {
    for (file in names(input_seurat[[run_v]])) {
      file_markers <- markers[[run_v]][[file]]
      file_markers_order <- ct_ordered[which(ct_ordered %in% unique(file_markers$cluster))]
      file_markers$cluster <- factor(file_markers$cluster,file_markers_order)
      file_markers <- file_markers[order(file_markers$cluster), ]
      hmp_top_order <- ct_ordered[which(ct_ordered %in% unique(file_markers$cluster))]
      Idents(input_seurat[[run_v]][[file]]) <- "ct"
      if (ord_levels=="yes") {
        levels(input_seurat[[run_v]][[file]]) <- hmp_top_order
      }
      hmp_top <- DoHeatmap(input_seurat[[run_v]][[file]], features = file_markers$gene, group.by = "ident", angle = 90, size = 3)
      pdf(paste0(plot_path, "hmp_all_", file, ".pdf"), height = 15)
      print(hmp_top  + labs(title = file))
      dev.off()
    }
  }
}

################################## Check SCENIC outputs in the SeuratObjects

# 11. Imports the SCENIC outputs, either from the GRNBoost2 or from the AUCell steps
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, which SCENIC output to import, 
    # and if we should follow the project order or not
  # Return: lists of SCENIC outputs

ImportSCENICresults <- function(main_dir, dis_type, res_folder, proj_order = "no") {
  if (dis_type!=F) {
    all_path <- paste0(main_dir, dis_type, "/", res_folder, "/sampled_100_cells_all/")
    projs <- c("GSE157827", "GSE174367", "PRJNA544731")
  } else {
    all_path <- paste0(main_dir, res_folder, "/sampled_100_cells_all/")
  }
  all2 <- list()
  runs <- c("_1", "_2", "_3")
  if (res_folder == "1_GRN") {
    all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
    all <- lapply(all_files,function(x) {
      read.table(file = x, 
                 sep = '\t', 
                 header = TRUE)
    })
    all <- lapply(1:length(all), function(x) all[[x]][order(-all[[x]]$importance),])
    names(all) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
    if (proj_order == "no") {
      for (run_v in runs) {
        only_1 <- names(all)[which(grepl(run_v, names(all)))]
        all_1 <- all[only_1]
        names(all_1) <- str_remove_all(names(all_1), ".tsv")
        all2 <- append(all2, list(all_1))
      }
      names(all2) <- runs
    } else if (proj_order == "yes") {
      for (proj_id in projs) {
        only_1 <- names(all)[which(grepl(proj_id, names(all)))]
        all_1 <- all[only_1]
        names(all_1) <- str_remove_all(names(all_1), ".tsv")
        all2 <- append(all2, list(all_1))
      }
      names(all2) <- projs
    }
  } else if (res_folder == "3_AUCell") {
    all_files <- list.files(path = all_path, pattern = "\\.csv$",full.names = T)
    all <- lapply(all_files,function(x) {
      read.table(file = x, 
                 sep = ',', 
                 header = TRUE)
    })
    names(all) <- list.files(path = all_path, pattern = "\\.csv$",full.names = F)
    if (proj_order == "yes") {
      for (run_v in runs) {
        only_1 <- names(all)[which(grepl(run_v, names(all)))]
        all_1 <- all[only_1]
        names(all_1) <- str_remove_all(names(all_1), ".csv")
        all2 <- append(all2, list(all_1))
      }
      names(all2) <- runs
    } else if (proj_order == "no") {
      for (proj_id in projs) {
        only_1 <- names(all)[which(grepl(proj_id, names(all)))]
        all_1 <- all[only_1]
        names(all_1) <- str_remove_all(names(all_1), ".csv")
        all2 <- append(all2, list(all_1))
      }
      names(all2) <- projs
    }
  }
  all2 <- all2[lapply(all2,length)>0]
  return(all2)
}

# 12. Generates a heatmap displaying the expression of either all the TFs/Tgs or only the top x in the SeuratObject built form the SCENIC inputs
  # Input: main directory where to save UMAP plots, if there is any sub-folder in the file tree, scenic GRN list, the SeuratObject, the order of the cts,
    # if we have to use a cutoff or not, and if we have to follow the order of the cts (sometimes the SeuratObject has issues with the ct order)
  # Return: nothing, saves the plot instead

TfTgSeuratExpression <- function(main_dir, dis_type, scenic_all, input_seurat, ct_ordered, cutoff = "no", ord_levels="yes") {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/Seurat/")
  } else {
    plot_path <- paste0(main_dir, "plots/Seurat/")
  }
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (run_v in names(scenic_all)) {
    for (file in names(scenic_all[[run_v]])) {
      if (cutoff == "no") {
        file_markers <- scenic_all[[run_v]][[file]]
      } else {
        file_markers <- scenic_all[[run_v]][[file]][1:cutoff, ]
      }
      hmp_top_order <- ct_ordered[which(ct_ordered %in% unique(input_seurat[[run_v]][[file]]@meta.data$ct))]
      Idents(input_seurat[[run_v]][[file]]) <- "ct"
      if (ord_levels=="yes") {
        levels(input_seurat[[run_v]][[file]]) <- hmp_top_order
      }
      for (k in c("TF", "target")) {
        hmp_top <- DoHeatmap(input_seurat[[run_v]][[file]], features = file_markers[, k], group.by = "ident", angle = 90, size = 3)
        if  (cutoff == "no") {
          pdf(paste0(plot_path, k,"_hmp_all_",  file, ".pdf"), height = 15)
          print(hmp_top  + labs(title = file))
          dev.off()
        } else {
          pdf(paste0(plot_path, k, "_hmp_top_", cutoff, "_", file, ".pdf"), height = 15)
          print(hmp_top  + labs(title = file))
          dev.off()
        }
      }
    }
  }
}

################################## Check TF and targets presence in GRNBoost2 outputs

# 13. Generates a df which contain information about how many TFs/Tgs are found in each run for each sex
  # Input: the GRN output, if there is any sub-folder in the file tree, which parametr to analyze (TF or targets), how many TF/targets to plot
  # Return: TF or target df

SCENICExtractGRN <- function(grn_out, dis_type, obj, cutoff) {
  grn_filt <- list()
  for (k in names(grn_out)) {
    grn_k <- list()
    for (i in names(grn_out[[k]])) {
      grn_k <- append(grn_k, list(grn_out[[k]][[i]][1:cutoff, ]))
    }
    names(grn_k) <- names(grn_out[[k]])
    grn_filt <- append(grn_filt, list(grn_k))
  }
  names(grn_filt) <- names(grn_out)
  grn_ids <- list()
  run_id <- vector()
  for (id_1 in names(grn_filt)) {
    for (id_2 in names(grn_filt[[id_1]])) {
      grns_v <- lapply(grn_filt[[id_1]], "[", c(obj))
      grn_ids <- append(grn_ids, list(grns_v[[id_2]][[obj]]))
      run_id <- c(run_id, id_2)
    }
  }
  names(grn_ids) <- run_id
  tot_grns <- unique(unlist(grn_ids))
  run_id <- rep(run_id, each=length(tot_grns))
  grn_id <- rep(tot_grns, length(run_id))
  grns <- data.frame(run_id, grn_id)
  grns$presence <- rep("no", nrow(grns))
  for (id in unique(grns$run_id)){
    grns[which(grns$run_id==id & grns$grn_id %in% grn_ids[[id]]),"presence"] <- "yes"
  }
  grns$presence <- factor(grns$presence, c("yes", "no"))
  grns <- grns[order(grns$presence), ]
  if (dis_type!=F) {
    cols_df <- c("proj", "sex", "run")
  } else {
    cols_df <- c("sex", "run")
  }
  grns <- separate(grns, run_id, into=cols_df, sep="_", remove=F)
  grns[cols_df] <- lapply(grns[cols_df], as.factor)
  grns <- grns[!duplicated(grns),]
  return(grns)
}

# 14. Generates and saves the plots for the TF/target expression in the runs and comapres the sexes
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree, scenic GRN list, if TF or targets are analyzed
  # Return: nothing, saves the plots instead

SCENICPlotGRN <- function(main_dir, dis_type, obj_grn, obj) {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/TF_TG/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    grnplot <- GRNPlot(obj_grn, obj, T)
  } else {
    plot_path <- paste0(main_dir, "plots/TF_TG/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    grnplot <- GRNPlot(obj_grn, obj)
  }
  pdf(paste0(plot_path, "hmp_", obj, ".pdf"))
  print(grnplot)
  dev.off()
}

# 15. Generates a heatmap displaying the rpesence of TF/targets in the three runs, F v M
  # Input: TF/target presnece df, TF or target to be analyzed, if the projects should be faceted
  # Return: presence heatmap

GRNPlot <- function(grn_data, obj, proj_facet=F) {
  grnplot <- ggplot(grn_data, aes(run, factor(grn_id, levels = rev(levels(factor(grn_id)))), fill=presence)) +
    geom_tile(color="#D3D3D3") +
    coord_fixed() +
    labs(x="Runs", y=obj, fill=paste0(obj, " found")) +
    {if (proj_facet) facet_wrap(~proj+sex, nrow = 1) } +
    {if (proj_facet==F) facet_wrap(~sex) } +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(grnplot)
}


################################## Plots the TFs and Targets expression in the original SeuratObject, split by celltypes and sex

# 16. Generates and saves to a CSV files which TF/targets are different between F and M
  # Input: the GRN output, if there is any sub-folder in the file tree, which parameter to analyze (TF or targets), how many TF/targets to plot
  # Return: nothing, saves the CSV to the output path

ExtractDiffGRN <- function(main_dir, dis_type, grns, obj) {
  if (dis_type!=F) {
    out_path <- paste0(main_dir, dis_type, "/5_outputs/")
    dir.create(out_path, showWarnings = F, recursive = T)
    grn_diff <- DISCOGRN(grns)
  } else {
    out_path <- paste0(main_dir, "5_outputs/")
    dir.create(out_path, showWarnings = F, recursive = T)
    grn_diff <- UCSCGRN(grns)
  }
  write.csv(grn_diff, paste0(out_path, "different_", obj, "_between_sexes.csv"))
}

# 17. Generates the different TF/target df for the DISCO datasets - multiple projects
  # Input: the GRN output
  # Return: the dataframe with the sex-biased TF/target

DISCOGRN <- function(grns) {
  grn_diff <- data.frame()
  for (proj_id in levels(grns$proj)) {
    gene_id <- vector()
    F_counts <- vector()
    M_counts <- vector()
    higher_expr <- vector()
    for (gene in unique(grns$grn_id)) {
      f_gene_count <- sum(grns[which(grns$proj==proj_id & grns$grn_id==gene & grns$sex == "F"), "presence"]=="yes")
      m_gene_count <- sum(grns[which(grns$proj==proj_id & grns$grn_id==gene & grns$sex == "M"), "presence"]=="yes")
      if (f_gene_count == 0 | m_gene_count == 0) {
        if (abs(f_gene_count - m_gene_count)>1) {
          gene_id <- c(gene_id , gene)
          F_counts<- c(F_counts, f_gene_count)
          M_counts<- c(M_counts, m_gene_count)
          if (f_gene_count - m_gene_count > 0 ) {
            higher_expr <- c(higher_expr, "F")
          } else {
            higher_expr <- c(higher_expr, "M")
          }
        }
      }
    }
    projs <- rep(proj_id, length(gene_id))
    grn_diff <- rbind(grn_diff, data.frame(projs, gene_id, F_counts, M_counts, higher_expr))
  }
  return(grn_diff)
}

# 18. Generates the different TF/target df for the UCSC dataset
  # Input: the GRN output
  # Return: the dataframe with the sex-biased TF/target

UCSCGRN <- function(grns) {
  gene_id <- vector()
  F_counts <- vector()
  M_counts <- vector()
  higher_expr <- vector()
  grns$sex <- str_replace_all(grns$sex, c("Female"="F", "Male"="M"))
  for (gene in unique(grns$grn_id)) {
    f_gene_count <- sum(grns[which(grns$grn_id==gene & grns$sex == "F"), "presence"]=="yes")
    m_gene_count <- sum(grns[which(grns$grn_id==gene & grns$sex == "M"), "presence"]=="yes")
    if (f_gene_count == 0 | m_gene_count == 0) {
      if (abs(f_gene_count - m_gene_count)>1) {
        gene_id <- c(gene_id , gene)
        F_counts<- c(F_counts, f_gene_count)
        M_counts<- c(M_counts, m_gene_count)
        if (f_gene_count - m_gene_count > 0 ) {
          higher_expr <- c(higher_expr, "F")
        } else {
          higher_expr <- c(higher_expr, "M")
        }
      }
    }
  }
  return(data.frame(gene_id, F_counts, M_counts, higher_expr))
}

# 19. Plot RdigePlots of sex-biased TFs and TGs
  # Input: main directory where to save the ridgeplots, the SeuratObject, the sex-biased TF or targets, 
    # how to group the cells in the SeuratObject (ct), and if TF or targets
  # Return:nothing, saves the plots instead

RidgeTFTG <- function(main_dir, seurat_obj, sex_biased_obj, obj_ident, obj) {
  if (length(sex_biased_obj)==0) {
    print(paste0("No unique ", obj, "s found"))
  } else {
    out_path <- paste0(main_dir, "plots/TF_TG/")
    dir.create(out_path, showWarnings = F, recursive = T)
    if (length(sex_biased_obj)==1) {
      pdf(paste0(out_path, obj, "_ridgeplot.pdf"))
      print(RidgePlot(seurat_obj, assay = "RNA", features = sex_biased_obj, group.by = obj_ident))
      dev.off()
    } else if (length(sex_biased_obj)>20) {
      sex_biased_obj_list <- split(sex_biased_obj, ceiling(seq_along(sex_biased_obj)/10))
      for (i in names(sex_biased_obj_list)) {
        pdf(paste0(out_path, obj, "_ridgeplot_", i, ".pdf"), width = 10)
        print(RidgePlot(seurat_obj, assay = "RNA", features = sex_biased_obj_list[[i]], group.by = obj_ident, stack = T))
        dev.off()
      }
    } else {
      pdf(paste0(out_path, obj, "_ridgeplot.pdf"), width = 10)
      print(RidgePlot(seurat_obj, assay = "RNA", features = sex_biased_obj, group.by = obj_ident, stack = T))
      dev.off()
    }
  }
}


################################## Check regulons presence in AUCell outputs

# 20. Genertaes a df with the presence of regulons across the runs and F v M
  # Input: the AUCell output, if there is any sub-folder in the file tree
  # Return: list of regulons presence dfs

SCENICExtractRegulons <- function(aucell_out, dis_type) {
  if (dis_type!=F) {
    regulons <- list()
    for (id in names(aucell_out)) {
      regulons_v <- lapply(aucell_out[[id]], "[", c("Regulon"))
      regulon_ids <- vector()
      run_id <- vector()
      for (id in names(regulons_v)){
        regulon_ids <- c(regulon_ids, regulons_v[[id]][["Regulon"]])
        run_id <- c(run_id, rep(id, length(regulons_v[[id]][["Regulon"]])))
      }
      tot_regs <- unique(regulon_ids)
      run_id <- rep(names(regulons_v), each=length(tot_regs))
      regulon_id <- rep(tot_regs, length(names(regulons_v)))
      reg_df <- data.frame(run_id, regulon_id)
      reg_df$presence <- rep("no", nrow(reg_df))
      for (id in names(regulons_v)){
        reg_df[which(reg_df$run_id==id & reg_df$regulon_id %in% regulons_v[[id]][["Regulon"]]),"presence"] <- "yes"
      }
      reg_df$presence <- factor(reg_df$presence, c("yes", "no"))
      reg_df <- reg_df[order(reg_df$presence), ]
      cols_df <- c("proj", "sex", "run")
      reg_df <- separate(reg_df, run_id, into=cols_df, sep="_", remove=F)
      reg_df[cols_df] <- lapply(reg_df[cols_df], as.factor)
      regulons_v <- list(reg_df)
      names(regulons_v) <- id
      regulons <- append(regulons, regulons_v)
    }
    names(regulons) <- names(aucell_out)
  } else {
    regulon_ids <- list()
    run_id <- vector()
    for (id_1 in names(aucell_out)) {
      for (id_2 in names(aucell_out[[id_1]])) {
        regulons_v <- lapply(aucell_out[[id_1]], "[", c("Regulon"))
        regulon_ids <- append(regulon_ids, list(regulons_v[[id_2]][["Regulon"]]))
        run_id <- c(run_id, id_2)
      }
    }
    names(regulon_ids) <- run_id
    tot_regs <- unique(unlist(regulon_ids))
    run_id <- rep(run_id, each=length(tot_regs))
    regulon_id <- rep(tot_regs, length(run_id))
    regulons <- data.frame(run_id, regulon_id)
    regulons$presence <- rep("no", nrow(regulons))
    for (id in unique(regulons$run_id)){
      regulons[which(regulons$run_id==id & regulons$regulon_id %in% regulon_ids[[id]]),"presence"] <- "yes"
    }
    regulons$presence <- factor(regulons$presence, c("yes", "no"))
    regulons <- regulons[order(regulons$presence), ]
    cols_df <- c("sex", "run")
    regulons <- separate(regulons, run_id, into=cols_df, sep="_", remove=F)
    regulons[cols_df] <- lapply(regulons[cols_df], as.factor)
  }
  return(regulons)
}

# 21. Generates and saves a heatmap displaying the presence of regulons across the runs and F v M
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree, list of regulons presence dfs
  # Return: nothing, saves the plot instead

SCENICPlotRegulons <- function(main_dir, dis_type, regulons) {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/Regulons/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    for (id in names(regulons)) {
      regplot <- RegulonsPlot(regulons[[id]], id)
      pdf(paste0(plot_path, id, "_hmp_regulons.pdf"))
      print(regplot)
      dev.off()
    }
  } else {
    plot_path <- paste0(main_dir, "plots/Regulons/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    regplot <- RegulonsPlot(regulons, "")
    pdf(paste0(plot_path, "hmp_regulons.pdf"))
    print(regplot)
    dev.off()
  }
}

# 22. Generates a heatmap displaying the presence of regulons across the runs and F v M
# Input: presence regulons df, id for the title (used if sub-groups are present)
# Return: plot

RegulonsPlot <- function(reg_data, id) {
  regplot <- ggplot(reg_data, aes(run, factor(regulon_id, levels = rev(levels(factor(regulon_id)))), fill=presence)) +
    geom_tile(color="#D3D3D3") +
    coord_fixed() +
    facet_wrap(~sex) +
    labs(x="Runs", y="Regulons", fill="Regulon found", title=id) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(regplot)
}

################################## Check output from individual cts and all - only on DISCO

# 23. Imports the GRN outputs when the celltypes were ran separately
# Input: main directory where to save the plots, if there is any sub-folder in the file tree
# Return: list of ct GRN outputs

SCENICct <- function(main_dir, dis_type) {
  ct_path <- paste0(main_dir, dis_type, "/1_GRN/sampled_100_cells/")
  ct_files <- list.files(path = ct_path, pattern = "\\.tsv$",full.names = T)
  cts <- lapply(ct_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(cts) <- list.files(path = ct_path, pattern = "\\.tsv$",full.names = F)
  names(cts) <- sapply(1:length(names(cts)), function(x) str_remove(names(cts)[x], ".tsv"))
  for (df in names(cts)) {
    cts[[df]]$TF_TG <- paste(cts[[df]]$TF, cts[[df]]$target, sep="_")
    cts[[df]]$TF_TG <- as.factor(cts[[df]]$TF_TG)
  }
  cts_sort <- lapply(1:length(names(cts)), function(x) cts[[x]][order(-cts[[x]]$importance),] )
  names(cts_sort) <- names(cts)
  return(cts_sort)
}

# 24. Imports the GRN outputs when the celltypes were ran together
# Input: main directory where to save the plots, if there is any sub-folder in the file tree
# Return: list of GRN outputs

SCENICall <- function(main_dir, dis_type) {
  all_path <- paste0(main_dir, dis_type, "/1_GRN/sampled_100_cells_all/")
  all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
  alls <- lapply(all_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(alls) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
  names(alls) <- sapply(1:length(names(alls)), function(x) str_remove(names(alls)[x], ".tsv"))
  for (df in names(alls)) {
    alls[[df]]$TF_TG <- paste(alls[[df]]$TF, alls[[df]]$target, sep="_")
    alls[[df]]$TF_TG <- as.factor(alls[[df]]$TF_TG)
  }
  alls_sort <- lapply(1:length(names(alls)), function(x) alls[[x]][order(-alls[[x]]$importance),] )
  names(alls_sort) <-names(alls) 
  return(alls_sort)
}

# 25. Imports both cts and all GRN outputs and combines them in a list
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree
  # Return: list of ct GRN and all GRN outputs combined in a list

SCENICoverlap <- function(main_dir, dis_type) {
  cts_sort <- SCENICct(main_dir, dis_type)
  alls_sort <- SCENICall(main_dir, dis_type)
  return(list(cts_sort, alls_sort))
}

# 26. Plots the % of overlap among runs, either in the all output or in each ct
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree, sorted list of both types of runs, the threshold of how many pairs to overlap,
    # the order of the celltypes
  # Return: nothing, saves the plot instead

PlotOverlapRuns <- function(main_dir, dis_type, sort_list, top_list, ct_list) {
  dir.create(paste0(main_dir, dis_type, "/plots/Overlap"), showWarnings = F, recursive = T)
  cts_sort <- sort_list[[1]]
  alls_sort <- sort_list[[2]]
  for (k in top_list) {
    print(k)
    if (k!="no") {
      k <- as.numeric(k)
      ct <- vector()
      tot <- vector()
      common <- vector()
      for (id in ct_list) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        ct <- c(ct, id)
        common <- c(common, length(intersect(sub_list[[1]][1:k, "TF_TG"], intersect(sub_list[[2]][1:k, "TF_TG"], sub_list[[3]][1:k, "TF_TG"]))))
      }
      tot <- rep(k, length(ct))
      overlap_ct <- data.frame(ct, tot, common)
      overlap_ct$perc <- overlap_ct$common*100 / overlap_ct$tot
      overlap_ct$ct <- as.factor(overlap_ct$ct)
      common_all <- length(intersect(alls_sort[[1]][1:k, "TF_TG"], intersect(alls_sort[[2]][1:k, "TF_TG"], alls_sort[[3]][1:k, "TF_TG"])))
      overlap_alls <- data.frame(c("All"), k, common_all, common_all*100 / as.numeric(k))
      colnames(overlap_alls) <- colnames(overlap_ct)
      overlap <- rbind(overlap_alls, overlap_ct)
      overlap$ct <- as.factor(overlap$ct)
      pdf(paste0(main_dir, dis_type, "/plots/Overlap/overlap_from_GRN_all_", k, "_TF_TG_pairs.pdf"))
      print(
        ggplot(overlap, aes(ct, perc, fill=ct)) +
          geom_bar(stat="identity") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    } else {
      ct <- vector()
      tot <- vector()
      common <- vector()
      for (id in ct_list) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        ct <- c(ct, id)
        tot <- c(tot, length(unique(sub_list[[1]][, "TF_TG"], sub_list[[2]][, "TF_TG"], sub_list[[3]][, "TF_TG"])))
        common <- c(common, length(intersect(sub_list[[1]][, "TF_TG"], intersect(sub_list[[2]][, "TF_TG"], sub_list[[3]][, "TF_TG"]))))
      }
      overlap_ct <- data.frame(ct, tot, common)
      overlap_ct$perc <- overlap_ct$common*100 / overlap_ct$tot
      overlap_ct$ct <- as.factor(overlap_ct$ct)
      tot_all <- length(unique(alls_sort[[1]][, "TF_TG"], alls_sort[[2]][, "TF_TG"], alls_sort[[3]][, "TF_TG"]))
      common_all <- length(intersect(alls_sort[[1]][, "TF_TG"], intersect(alls_sort[[2]][, "TF_TG"], alls_sort[[3]][, "TF_TG"])))
      overlap_alls <- data.frame(c("All"), tot_all, common_all, common_all*100 / tot_all)
      colnames(overlap_alls) <- colnames(overlap_ct)
      overlap <- rbind(overlap_alls, overlap_ct)
      overlap$ct <- as.factor(overlap$ct)
      pdf(paste0(main_dir, dis_type, "/plots/Overlap/unsorted_overall_from_GRN_all_TF_TG_pairs.pdf"))
      print(
        ggplot(overlap, aes(ct, perc, fill=ct)) +
          geom_bar(stat="identity") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    }
  }
}

################################## Check overlap with scGRNom

# 27. Plots the overlap between the single ct runs and the results from scGRNom
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree, sorted list of both types of runs, the pairs from scGRNom,
    # the common annotation vector, sorted comparison
  # Return: nothing, saves plot instead

PlotscGRNomOverlap <- function(main_dir, dis_type, sort_list, top_scGRNom, ct_scGRNom, flag) {
  dir.create(paste0(main_dir, dis_type, "/plots/scGRNom"), showWarnings = F, recursive = T)
  cts_sort <- sort_list[[1]]
  alls_sort <- sort_list[[2]]
  for (f in flag) {
    if (f == "no") {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (id in names(ct_scGRNom)) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        names(sub_list) <- sapply(1:length(names(sub_list)), function(x) str_replace(names(sub_list)[x], id, ct_scGRNom[[id]]))
        sub_scGRNom_names <- names(scGRNom)[which(grepl(ct_scGRNom[[id]], names(scGRNom)))]
        sub_scGRNom <- scGRNom[sub_scGRNom_names]
        for (scGRNom_id in names(sub_scGRNom)) {
          for (id_sub in names(sub_list)) {
            sort_sub <- sub_list[[id_sub]][order(-sub_list[[id_sub]]$importance),]
            sort_sub <- sort_sub[1:top_scGRNom[[scGRNom_id]], ]
            scGRNom_groups <- c(scGRNom_groups, paste(id_sub, scGRNom_id, sep="_"))
            tot_scGRNom <- c(tot_scGRNom, length(unique(sort_sub$TF_TG, sub_scGRNom[scGRNom_id]$TF_TG)))
            common_scGRNom <- c(common_scGRNom, length(intersect(sort_sub$TF_TG, sub_scGRNom[[scGRNom_id]]$TF_TG)))
          }
        }
      }
      overlap_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_scGRNom$perc <- overlap_scGRNom$common_scGRNom * 100 / overlap_scGRNom$tot_scGRNom
      overlap_scGRNom$scGRNom_groups <- str_replace_all(overlap_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_scGRNom <- separate(overlap_scGRNom, scGRNom_groups, into=c("disco", "run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("disco", "run", "GRN_ct", "openchrom")
      overlap_scGRNom[col_factors] <- lapply(overlap_scGRNom[col_factors], as.factor)
      pdf(paste0(main_dir, dis_type, "/plots/scGRNom/unsorted_overall_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black")))
      dev.off()
    } else {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (id in names(ct_scGRNom)) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        names(sub_list) <- sapply(1:length(names(sub_list)), function(x) str_replace(names(sub_list)[x], id, ct_scGRNom[[id]]))
        sub_scGRNom_names <- names(scGRNom)[which(grepl(ct_scGRNom[[id]], names(scGRNom)))]
        sub_scGRNom <- scGRNom[sub_scGRNom_names]
        for (scGRNom_id in names(sub_scGRNom)) {
          for (id_sub in names(sub_list)) {
            sort_sub <- sub_list[[id_sub]][order(-sub_list[[id_sub]]$importance),]
            sort_sub <- sort_sub[1:top_scGRNom[[scGRNom_id]], ]
            scGRNom_groups <- c(scGRNom_groups, paste(id_sub, scGRNom_id, sep="_"))
            tot_scGRNom <- c(tot_scGRNom, top_scGRNom[[scGRNom_id]])
            common_scGRNom <- c(common_scGRNom, length(intersect(sort_sub$TF_TG, sub_scGRNom[[scGRNom_id]]$TF_TG)))
          }
        }
      }
      overlap_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_scGRNom$perc <- overlap_scGRNom$common_scGRNom * 100 / overlap_scGRNom$tot_scGRNom
      overlap_scGRNom$scGRNom_groups <- str_replace_all(overlap_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_scGRNom <- separate(overlap_scGRNom, scGRNom_groups, into=c("disco", "run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("disco", "run", "GRN_ct", "openchrom")
      overlap_scGRNom[col_factors] <- lapply(overlap_scGRNom[col_factors], as.factor)
      pdf(paste0(main_dir, dis_type, "/plots/scGRNom/unsorted_overlap_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black")))
      dev.off()
    }
  }
  for (f in flag) {
    if (f=="no") {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (a in names(alls_sort)) {
        for (c in names(scGRNom)) {
          scGRNom_groups <- c(scGRNom_groups, paste(a, c, sep="_"))
          a_sort <- alls_sort[[a]][1:top_scGRNom[[c]], "TF_TG"]
          common_scGRNom <- c(common_scGRNom, length(intersect(a_sort, sub_scGRNom[[c]]$TF_TG)))
          tot_scGRNom <- c(tot_scGRNom, length(unique(a_sort, sub_scGRNom[[c]]$TF_TG)))
        }
      }
      overlap_alls_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_alls_scGRNom$scGRNom_groups <- str_replace_all(overlap_alls_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_alls_scGRNom <- separate(overlap_alls_scGRNom, scGRNom_groups, into=c("run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("run", "GRN_ct", "openchrom")
      overlap_alls_scGRNom[col_factors] <- lapply(overlap_alls_scGRNom[col_factors], as.factor)
      overlap_alls_scGRNom$perc <- overlap_alls_scGRNom$common_scGRNom * 100 / overlap_alls_scGRNom$tot_scGRNom
      pdf(paste0(main_dir, dis_type, "/plots/scGRNom/unsorted_overall_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    } else {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (a in names(alls_sort)) {
        for (c in names(scGRNom)) {
          scGRNom_groups <- c(scGRNom_groups, paste(a, c, sep="_"))
          tot_scGRNom <- c(tot_scGRNom, top_scGRNom[[c]])
          a_sort <- alls_sort[[a]][1:top_scGRNom[[c]], "TF_TG"]
          common_scGRNom <- c(common_scGRNom, length(intersect(a_sort, sub_scGRNom[[c]]$TF_TG)))
        }
      }
      overlap_alls_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_alls_scGRNom$scGRNom_groups <- str_replace_all(overlap_alls_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_alls_scGRNom <- separate(overlap_alls_scGRNom, scGRNom_groups, into=c("run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("run", "GRN_ct", "openchrom")
      overlap_alls_scGRNom[col_factors] <- lapply(overlap_alls_scGRNom[col_factors], as.factor)
      overlap_alls_scGRNom$perc <- overlap_alls_scGRNom$common_scGRNom * 100 / overlap_alls_scGRNom$tot_scGRNom
      pdf(paste0(main_dir, dis_type, "/plots/scGRNom/unsorted_overlap_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    }
  }
}

################################## Check number of TF-TG pairs between M and F within the same project

# 28. Creates the TF-target pairs from the GRN outputs
  # Input: grn output, if there is any sub-folder in the file tree
  # Return: grn output with added column in all element of list

GRNpairs <- function(all_grn, dis_type) {
  if (dis_type!=F) {
    for (proj_id in names(all_grn)) {
      for (run in names(all_grn[[proj_id]])) {
        all_grn[[proj_id]][[run]]$pairs <- paste(all_grn[[proj_id]][[run]]$TF,  all_grn[[proj_id]][[run]]$target, sep="_")
      }
    }
  } else {
    for (run in names(all_grn)) {
      for (sex in names(all_grn[[run]])) {
        all_grn[[run]][[sex]]$pairs <- paste(all_grn[[run]][[sex]]$TF,  all_grn[[run]][[sex]]$target, sep="_")
      }
    }
  }
  return(all_grn)
}

# 29. Calculates overlap of TF-target pairs between one run of one sex and the three runs of the other sex
  # Input: GRN output, if there is any sub-folder in the file tree, a cutoff for how many pairs to compare, and if Velmehsev, change the sex values
  # Return: list of pair count dfs

SCENICOverlapTfTg <-  function(all_grn, dis_type=F, threshold="no", analysis_type="no") {
  all_grn <- GRNpairs(all_grn, dis_type)
  id <- vector()
  overlap <- vector()
  count <- vector()
  if (dis_type!=F) {
    sexes <- c("_F_", "_M_")
  } else {
    if (analysis_type=="Velmeshev") {
      sexes <- c("Female_", "Male_")
    } else {
      sexes <- c("F_", "M_")
    }
    all_grn <- unlist(all_grn, recursive = F)
  }
  for (sex in sexes) {
    if (dis_type!=F) {
      for (proj_id in names(all_grn)) {
        only_1 <- names(all_grn[[proj_id]])[which(grepl(sex, names(all_grn[[proj_id]])))]
        sex_proj <- all_grn[[proj_id]][only_1]
        other_sex <- sexes[which(sexes!=sex)]
        other_sex_proj <- names(all_grn[[proj_id]])[which(grepl(other_sex, names(all_grn[[proj_id]])))]
        other_sex_proj <- all_grn[[proj_id]][other_sex_proj]
        if (threshold=="no") {
          common_pairs <- sapply(sex_proj, "[", c("pairs"))
          common_pairs <- unique(unlist(common_pairs))
        } else {
          sex_proj <- do.call(rbind, sex_proj)
          common_pairs <- sex_proj[order(-sex_proj$importance),"pairs"]
          common_pairs <- unique(common_pairs)[1:threshold]
        }
        sex_pairs <- sapply(other_sex_proj, "[[", c("pairs"))
        sex_pairs[[sex]] <- common_pairs
        sex_pairs <- as.data.frame(t(table(unlist(sex_pairs))))
        sex_pairs[,1] <- NULL
        colnames(sex_pairs) <- c("pairs", "count")
        for (k in 1:length(unique(sex_pairs$count))) {
          id <- c(id, paste0(proj_id, stri_replace_last_fixed(sex, "_", "")))
          overlap <- c(overlap, k-1)
          count <- c(count, sum(sex_pairs$count==k))
        }
      }  
    } else {
      only_1 <- names(all_grn)[which(grepl(sex, names(all_grn)))]
      sex_df <- all_grn[only_1]
      other_sex <- sexes[which(sexes!=sex)]
      other_sex_df <- names(all_grn)[which(grepl(other_sex, names(all_grn)))]
      other_sex_df <- all_grn[other_sex_df]
      if (threshold=="no") {
      common_pairs <- sapply(sex_df, "[", c("pairs"))
      common_pairs <- unique(unlist(common_pairs))
      } else {
        sex_df <- do.call(rbind, sex_df)
        common_pairs <- sex_df[order(-sex_df$importance),"pairs"]
        common_pairs <- unique(common_pairs)[1:threshold]
      }
      sex_pairs <- sapply(other_sex_df, "[", c("pairs"))
      sex_pairs[[sex]] <- common_pairs 
      sex_pairs <- as.data.frame(t(table(unlist(sex_pairs))))
      sex_pairs[,1] <- NULL
      colnames(sex_pairs) <- c("pairs", "count")
      for (k in 1:length(unique(sex_pairs$count))) {
        id <- c(id, str_replace_all(sex, "_", ""))
        overlap <- c(overlap, k-1)
        count <- c(count, sum(sex_pairs$count==k))
      }
    }
    tot_k <- length(unique(sex_pairs$count)) - 1
  }
  df_counts <- data.frame(id, overlap, count)
  df_counts[which(df_counts$overlap==0), "overlap"] <- "None"
  if (dis_type!=F) {
    df_counts <- separate(df_counts, id, into=c("proj", "sex"), remove = F, sep="_")
    df_counts[c("id", "proj", "sex", "overlap")] <- lapply(df_counts[c("id", "proj", "sex", "overlap")], as.factor)
  } else {
    df_counts[c("id", "overlap")] <- lapply(df_counts[c("id", "overlap")], as.factor)
  }
  df_counts$overlap <- factor(df_counts$overlap, c("None", seq(1,tot_k)))
  return(df_counts)
}

# 30. Plots the number of overlapping TF/target pairs both as bar lot and %
  # Input: main directory where to save the plots, if there is any sub-folder in the file tree, the overlap counts df, if any threshold was used
  # Return: nothing, saves the plot instead

SCENICPlotOverlapTfTg <- function(main_dir, dis_type, df_counts, threshold = F) {
  if (dis_type!=F) {
    plot_path <- paste0(main_dir, dis_type, "/plots/TF_TG/")
  } else {
    plot_path <- paste0(main_dir, "plots/TF_TG/")
  }
  dir.create(plot_path, showWarnings = F, recursive = T)
  p_dodge <- ggplot(df_counts, aes(id, count, fill=overlap)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(x="Reference TF-TG pairs", y="Counts of TF-TG pairs", fill="Overlap") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  p_fill <- ggplot(df_counts, aes(id, count, fill=overlap)) +
    geom_bar(stat = "identity", position = "fill") +
    labs(x="Reference TF-TG pairs", y="% of TF-TG pairs", fill="Overlap") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  TfTgplots <- ggarrange(p_dodge, p_fill, common.legend = T, legend = "bottom")
  if (threshold!=F) {
    pdf(paste0(plot_path, "pairs_overlap_among sexes_thresh_", threshold, ".pdf"))
    print(TfTgplots)
    dev.off()
  } else {
    pdf(paste0(plot_path, "pairs_overlap_among sexes.pdf"))
    print(TfTgplots)
    dev.off()
  }
}
