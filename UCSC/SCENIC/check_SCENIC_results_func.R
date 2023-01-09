library(ggplot2)
library(readxl)
library(stringr)
library(stringi)
library(tidyr)
library(reshape2)
library(Seurat)
library(SeuratObject)
library(ggpubr)
library(S4Vectors)
library(purrr)
library(dplyr)
library(Seurat)
library(matrixStats)

################################## Check randomly sampled SCENIC inputs

SCENICInputSeurat <- function(main_dir, dis_type, run_v) {
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
  #names(input_dfs) <- str_remove_all(names(input_dfs), "_1")
  for (id in names(input_dfs)) {
    #colnames(input_dfs[[id]]) <- str_replace_all(colnames(input_dfs[[id]]), "[.]", "-")
    rownames(input_dfs[[id]]) <- input_dfs[[id]]$Genes
    input_dfs[[id]]$Genes <- NULL
  #  sub_info <- subset(cell_info_df, cell_id %in% colnames(input_dfs[[id]])[-1])
  #  ct_info <- sub_info$ct
  #  colnames(input_dfs[[id]]) <- c("Genes", ct_info)
  }
  return(input_dfs)
}

# runs the Seurat clustering from expression matrix (follows tutorial)
SCENICSeuratPlots <- function(expr_matrix_list, metadata_id, id) {
  id_seurat <- CreateSeuratObject(as.matrix(expr_matrix_list[[id]]), project = id,
                                    assay = "RNA", meta.data = metadata_id)
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  id_seurat[["percent.mt"]] <- PercentageFeatureSet(id_seurat, pattern = "^MT-")
  # Visualize QC metrics as a violin plot
  #VlnPlot(id_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  #plot1 <- FeatureScatter(id_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
  #plot2 <- FeatureScatter(id_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #plot1 + plot2
  #id_seurat <- subset(id_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  id_seurat <- subset(id_seurat, subset = percent.mt < 5)
  id_seurat <- NormalizeData(id_seurat)
  id_seurat <- FindVariableFeatures(id_seurat, selection.method = "vst", nfeatures = 2000)
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(id_seurat), 10)
  # plot variable features with and without labels
  #plot1 <- VariableFeaturePlot(id_seurat)
  #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #plot1 + plot2
  all.genes <- rownames(id_seurat)
  id_seurat <- ScaleData(id_seurat, features = all.genes)
  id_seurat <- RunPCA(id_seurat, features = VariableFeatures(object = id_seurat))
  # Examine and visualize PCA results a few different ways
  print(id_seurat[["pca"]], dims = 1:5, nfeatures = 5)
  VizDimLoadings(id_seurat, dims = 1:2, reduction = "pca")
  #DimPlot(id_seurat, reduction = "pca") 
  #DimHeatmap(id_seurat, dims = 1, cells = 500, balanced = TRUE)
  p1 <- DimHeatmap(id_seurat, dims = 1:15, cells = 500, balanced = TRUE)
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More
  # approximate techniques such as those implemented in ElbowPlot() can be used to reduce
  # computation time
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

# after deciding the number of PCAs to use, plots the final UMAP
SCENICClustering <- function(main_dir, dis_type, id_seurat, cluster_num, ct_ordered, plot_flag = "no", file_out_name) {
  id_seurat <- FindNeighbors(id_seurat, dims = 1:cluster_num)
  id_seurat <- FindClusters(id_seurat, resolution = 0.5)
  # Look at cluster IDs of the first 5 cells
  #head(Idents(id_seurat), 5)
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
  # 'umap-learn')
  id_seurat <- RunUMAP(id_seurat, dims = 1:cluster_num)
  # note that you can set `label = TRUE` or use the LabelClusters function to help label
  # individual clusters
  if (plot_flag == "yes") {
    SCENICUmap(main_dir, dis_type, id_seurat, ct_ordered, file_out_name)
  }
  return(id_seurat)
}

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


Filter_gene <- function(order.gene.df, pval, FC) {
  logFC <- log2(1/FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  return(gene.sig)
}

HmpSCENIC <- function(main_dir, dis_type, input_seurat, top10, ct_ordered, ord_level="yes") {
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
      if (ord_level=="yes") {
        levels(input_seurat[[run_v]][[file]]) <- hmp_top_order
      }
      hmp_top <- DoHeatmap(input_seurat[[run_v]][[file]], features = topgenes$genes, group.by = "ident", angle = 90, size = 3)
      pdf(paste0(plot_path, "hmp_top_10_", file, ".pdf"), height = 15)
      print(hmp_top  + labs(title = file))
      dev.off()
    }
  }
}

HmpSCENICAll <- function(main_dir, dis_type, input_seurat, markers, ct_ordered, ord_level="yes") {
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
      if (ord_level=="yes") {
        levels(input_seurat[[run_v]][[file]]) <- hmp_top_order
      }
      hmp_top <- DoHeatmap(input_seurat[[run_v]][[file]], features = file_markers$gene, group.by = "ident", angle = 90, size = 3)
      pdf(paste0(plot_path, "hmp_all_", file, ".pdf"), height = 15)
      print(hmp_top  + labs(title = file))
      dev.off()
    }
  }
  #return(hmp_top)
}

SCENICresultsSeurat <- function(main_dir, dis_type, res_folder, proj_order = "no") {
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

SCENICTfTg <- function(main_dir, dis_type, scenic_all, input_seurat, ct_ordered, cutoff = "no", ord_levels="yes") {
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

# Plot RdigePlots of unique TFs and TGs
RidgeTFTG <- function(main_dir, seurat_obj, shared_features, obj_ident, obj) {
  if (length(shared_features)==0) {
    print(paste0("No unique ", obj, "s found"))
  } else {
    out_path <- paste0(main_dir, "plots/TF_TG/")
    dir.create(out_path, showWarnings = F, recursive = T)
    if (length(shared_features)==1) {
      pdf(paste0(out_path, obj, "_ridgeplot.pdf"))
      print(RidgePlot(seurat_obj, assay = "RNA", features = shared_features, group.by = obj_ident))
      dev.off()
    } else if (length(shared_features)>20) {
      shared_features_list <- split(shared_features, ceiling(seq_along(shared_features)/20))
      for (i in names(shared_features_list)) {
        pdf(paste0(out_path, obj, "_ridgeplot_", i, ".pdf"), width = 10)
        print(RidgePlot(seurat_obj, assay = "RNA", features = shared_features_list[[i]], group.by = obj_ident, stack = T))
        dev.off()
      }
    } else {
      pdf(paste0(out_path, obj, "_ridgeplot.pdf"), width = 10)
      print(RidgePlot(seurat_obj, assay = "RNA", features = shared_features, group.by = obj_ident, stack = T))
      dev.off()
    }
  }
}



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

################################## Check regulons presence in AUCell outputs

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

################################## Check output from individual cts and all

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

SCENICoverlap <- function(main_dir, dis_type) {
  cts_sort <- SCENICct(main_dir, dis_type)
  alls_sort <- SCENICall(main_dir, dis_type)
  return(list(cts_sort, alls_sort))
}

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

################################## Check expression of TFs and TGs in original data


#SCENICresults <- function(main_dir, dis_type) {
#  dir.create(paste0(main_dir, dis_type, "/plots"), showWarnings = F, recursive = T)
#  all_path <- paste0(main_dir, dis_type, "/1_GRN/sampled_100_cells_all/")
#  all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
#  all <- lapply(all_files,function(x) {
#    read.table(file = x, 
#               sep = '\t', 
#               header = TRUE)
#  })
#  names(all) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
#  only_1 <- names(all)[which(grepl("_1", names(all)))]
#  all_1 <- all[only_1]
#  names(all_1) <- str_remove_all(names(all_1), "_1.tsv")
#  all_sort <- lapply(1:length(names(all_1)), function(x) all_1[[x]][order(-all_1[[x]]$importance),])
#  names(all_sort) <- names(all_1)
#  return(all_sort)
#}

SCENICInput <- function(main_dir, dis_type, cell_info_df) {
  input_dfs_path <-  paste0(main_dir, dis_type, "/0_input_dfs/sampled_100_cells_all")
  input_dfs_files <- list.files(path = input_dfs_path, full.names = T)
  input_dfs <- lapply(input_dfs_files,function(x) {
    read.table(file = x, 
               sep = ',', 
               header = TRUE)
  })
  names(input_dfs) <- list.files(path = input_dfs_path, full.names = F)
  names(input_dfs) <- str_remove_all(names(input_dfs), ".csv")
  only_1 <- names(input_dfs)[which(grepl("_1", names(input_dfs)))]
  input_dfs <- input_dfs[only_1]
  names(input_dfs) <- str_remove_all(names(input_dfs), "_1")
  for (id in names(input_dfs)) {
    colnames(input_dfs[[id]]) <- str_replace_all(colnames(input_dfs[[id]]), "[.]", "-")
    sub_info <- subset(cell_info_df, cell_id %in% colnames(input_dfs[[id]])[-1])
    ct_info <- sub_info$ct
    colnames(input_dfs[[id]]) <- c("Genes", ct_info)
  }
  return(input_dfs)
}

PlotTfTg <- function(main_dir, dis_type, all_sort, input_dfs, ct_ordered, top_TF){
  dir.create(paste0(main_dir, dis_type, "/plots/TF_TG"), showWarnings = F, recursive = T)
  for (id in names(all_sort)) {
    for (k in top_TF) {
      tf_id <- all_sort[[id]][1:k, "TF"]
      sub_tf <- subset(input_dfs[[id]], Genes %in% tf_id)
      sub_tf <-  melt(sub_tf, id.vars = "Genes")
      colnames(sub_tf) <- c("TF", "ct", "value")
      sub_tf$ct <- as.character(sub_tf$ct)
      sub_tf$ct <- str_remove_all(sub_tf$ct, "[.][:digit:]+")
      sub_tf[c("TF", "ct")] <- lapply(sub_tf[c("TF", "ct")], as.factor)
      sub_tf_order <- ct_ordered[which(ct_ordered %in% levels(sub_tf$ct))]
      sub_tf$ct <- factor(sub_tf$ct, sub_tf_order)
      sub_tf <- sub_tf[order(sub_tf$ct), ]
      sub_tf$log2_value <- log2(sub_tf$value)
      sub_tf[which(sub_tf$log2_value=="-Inf"),"log2_value"] <- 0
      pdf(paste0(main_dir, dis_type, "/plots/TF_TG/TFs_hmp_",k, "_top_", id, ".pdf"))
      print(
        #ggplot(sub_tf, aes(ct, TF, fill=log2_value, color="")) +
        ggplot(sub_tf, aes(ct, TF, fill=log2_value)) +
          geom_tile(color="#D3D3D3") + 
          coord_fixed() +
          scale_fill_gradient2(low="blue", high="red", mid = "white", na.value = "#D3D3D3") +
          #scale_colour_manual(values=NA) + 
          #guides(colour=guide_legend("NA", override.aes=list(fill="#D3D3D3"))) +
          labs(y="TFs", fill="Log2 RNA Expression", title=id) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "right", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
      tg_id <- all_sort[[id]][1:k, "target"]
      sub_tg <- subset(input_dfs[[id]], Genes %in% tg_id)
      if (nrow(sub_tg) > 0) {
        sub_tg <- melt(sub_tg, id.vars = "Genes")
        colnames(sub_tg) <- c("TG", "ct", "value")
        sub_tg$ct <- as.character(sub_tg$ct)
        sub_tg$ct <- str_remove_all(sub_tg$ct, "[.][:digit:]+")
        sub_tg[c("TG", "ct")] <- lapply(sub_tg[c("TG", "ct")], as.factor)
        sub_tg_order <- ct_ordered[which(ct_ordered %in% levels(sub_tg$ct))]
        sub_tg$ct <- factor(sub_tg$ct, sub_tg_order)
        sub_tg <- sub_tg[order(sub_tg$ct), ]
        sub_tg$log2_value <- log2(sub_tg$value)
        sub_tg[which(sub_tg$log2_value=="-Inf"),"log2_value"] <- 0
        if (k==100) {
          pdf(paste0(main_dir, dis_type, "/plots/TF_TG/TGs_hmp_",k, "_top_", id, ".pdf"),
              height=11)
        } else {
          pdf(paste0(main_dir, dis_type, "/plots/TF_TG/TGs_hmp_",k, "_top_", id, ".pdf"),
              height=15)
        }
        print(
          #ggplot(sub_tg, aes(ct, TG, fill=log2_value, colour="")) +
          ggplot(sub_tg, aes(ct, TG, fill=log2_value)) +
            geom_tile(color="#D3D3D3") + 
            coord_fixed() +
            scale_fill_gradient2(low="blue", high="red", mid = "white", na.value = "#D3D3D3") +
            #scale_colour_manual(values=NA) + 
            #guides(colour=guide_legend("NA", override.aes=list(fill="#D3D3D3"))) +
            labs(y="TGs", fill="Log2 RNA Expression", title=id) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  plot.title = element_text(size=12, face="bold", colour = "black"),
                  axis.line = element_line(colour = "black"),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                  axis.ticks.x=element_blank(),
                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                  legend.position = "right", 
                  legend.title = element_text(size=12, face="bold", colour = "black"))
        )
        dev.off()
      }
    }
  }
}


################################## Check number of TF-TG pairs between M and F within the same project

SCENICAddTFTG <- function(all_grn, dis_type) {
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

SCENICOverlapTfTg <-  function(all_grn, dis_type, analysis_type="no") {
  all_grn <- SCENICAddTFTG(all_grn, dis_type)
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
        common_pairs <- sapply(sex_proj, "[", c("pairs"))
        common_pairs <- unique(unlist(common_pairs))
        other_sex <- sexes[which(sexes!=sex)]
        other_sex_proj <- names(all_grn[[proj_id]])[which(grepl(other_sex, names(all_grn[[proj_id]])))]
        other_sex_proj <- all_grn[[proj_id]][other_sex_proj]
        sex_pairs <- sapply(other_sex_proj, "[[", c("pairs"))
        sex_pairs[[sex]] <- common_pairs 
        sex_pairs <- as.data.frame(t(table(unlist(sex_pairs))))
        sex_pairs[,1] <- NULL
        colnames(sex_pairs) <- c("pairs", "count")
        for (k in 1:length(unique(sex_pairs$count))) {
          id <- c(id, paste0(proj_id,stri_replace_last_fixed(sex, "_", "")))
          overlap <- c(overlap, k-1)
          count <- c(count, sum(sex_pairs$count==k))
        }
      }  
      tot_k <- length(unique(sex_pairs$count)) - 1
    } else {
      only_1 <- names(all_grn)[which(grepl(sex, names(all_grn)))]
      sex_df <- all_grn[only_1]
      common_pairs <- sapply(sex_df, "[", c("pairs"))
      common_pairs <- unique(unlist(common_pairs))
      other_sex <- sexes[which(sexes!=sex)]
      other_sex_df <- names(all_grn)[which(grepl(other_sex, names(all_grn)))]
      other_sex_df <- all_grn[other_sex_df]
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
    tot_k <- length(unique(sex_pairs$count)) - 1
    }
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


SCENICPlotOverlapTfTg <- function(main_dir, dis_type, df_counts) {
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
  pdf(paste0(plot_path, "pairs_overlap_among sexes.pdf"))
  print(TfTgplots)
  dev.off()
}


###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################


################ OLDER SCRIPTS AND PLOTS

########## Plot cell types back to umap/tsne

##### Solution 1 - create a Seurat object for each expression matrix

# Normal

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

#all_norm_final <- readRDS(paste0(main, sub_disease[3], "/seurat_files.rds"))

##### Solution 2 - subset from original disco 

#meta_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[1]])),]
#rownames(meta_id) <- meta_id$cell_id
#meta_id$cell_id <- NULL


#disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

#disco_filt$cell_id <- rownames(disco_filt@meta.data)
#disco_filt@meta.data$ct <- str_replace_all(disco_filt@meta.data$ct, "/", "_")

#for (norm_id in names(norm_input_seurat)) {
#  metadata_id <- cell_info[which(cell_info$cell_id %in% colnames(norm_input_seurat[[norm_id]])),]
#  rownames(metadata_id) <- metadata_id$cell_id
#  metadata_id$cell_id <- NULL
#  id_sub <- subset(disco_filt, cell_id %in% rownames(metadata_id))
#  umap_id_sub <- DimPlot(id_sub, reduction = "umap", group.by = "ct")
#  umap_order <- ct_order[which(ct_order %in% levels(umap_id_sub$data$ct))]
#  umap_id_sub$data$ct <- factor(umap_id_sub$data$ct, umap_order)
#  umap_id_sub$data <- umap_id_sub$data[order(umap_id_sub$data$ct), ]
#  pdf(paste0(main, sub_disease[3], "/3_plots/UMAP_", norm_id, "_disco_subset.pdf"))
#  print(umap_id_sub  + labs(title = norm_id))
#  dev.off()
#}


########## Heatmap expression of SCENIC TFs and TGs 

#top_TFs <- c(100, 200)
#norm_sort <- SCENICresultsSeurat(main, sub_disease[3], "1_GRN")
#norm_input <- SCENICInput(main, sub_disease[3], cell_info)
#PlotTfTg(main, sub_disease[3], norm_sort, norm_input, ct_order, top_TFs)