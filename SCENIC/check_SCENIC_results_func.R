library(ggplot2)
library(readxl)
library(stringr)
library(tidyr)
library(reshape2)
library(Seurat)
library(SeuratObject)
library(ggpubr)
library(S4Vectors)

SCENICInputSeurat <- function(main_dir, dis_type, run_v) {
  input_dfs_path <-  paste0(main_dir, dis_type, "/0_input_dfs/sampled_100_cells_all")
  input_dfs_files <- list.files(path = input_dfs_path, full.names = T)
  input_dfs <- lapply(input_dfs_files,function(x) {
    read.table(file = x, 
               sep = ',', 
               header = TRUE)
  })
  names(input_dfs) <- list.files(path = input_dfs_path, full.names = F)
  names(input_dfs) <- str_remove_all(names(input_dfs), ".csv")
  only_1 <- names(input_dfs)[which(grepl(run_v, names(input_dfs)))]
  input_dfs <- input_dfs[only_1]
  #names(input_dfs) <- str_remove_all(names(input_dfs), "_1")
  for (id in names(input_dfs)) {
    colnames(input_dfs[[id]]) <- str_replace_all(colnames(input_dfs[[id]]), "[.]", "-")
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
SCENICUmap <- function(main_dir, dis_type, id_seurat, cluster_num, ct_ordered) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  id_seurat <- FindNeighbors(id_seurat, dims = 1:cluster_num)
  id_seurat <- FindClusters(id_seurat, resolution = 0.5)
  # Look at cluster IDs of the first 5 cells
  #head(Idents(id_seurat), 5)
  # If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
  # 'umap-learn')
  id_seurat <- RunUMAP(id_seurat, dims = 1:cluster_num)
  # note that you can set `label = TRUE` or use the LabelClusters function to help label
  # individual clusters
  umap_id <- DimPlot(id_seurat, reduction = "umap", group.by = "ct")
  umap_order <- ct_ordered[which(ct_ordered %in% levels(umap_id$data$ct))]
  umap_id$data$ct <- factor(umap_id$data$ct, umap_order)
  umap_id$data <- umap_id$data[order(umap_id$data$ct), ]
  pdf(paste0(main_dir, dis_type, "/3_plots/UMAP_", unique(id_seurat$orig.ident), ".pdf"))
  print(umap_id  + labs(title = unique(id_seurat$orig.ident)))
  dev.off()
}



SCENICct <- function(main_dir, dis_type) {
  ct_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells/")
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
  all_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells_all/")
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
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  cts_sort <- sort_list[[1]]
  alls_sort <- sort_list[[2]]
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
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
      pdf(paste0(main_dir, dis_type, "/3_plots/overlap_from_GRN_all_", k, "_TF_TG_pairs.pdf"))
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
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_all_TF_TG_pairs.pdf"))
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

PlotscGRNomOverlap <- function(main_dir, dis_type, sort_list, top_scGRNom, ct_scGRNom, flag) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
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
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
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
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overlap_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
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
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
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
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overlap_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
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

SCENICresults <- function(main_dir, dis_type) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  all_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells_all/")
  all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
  all <- lapply(all_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(all) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
  only_1 <- names(all)[which(grepl("_1", names(all)))]
  all_1 <- all[only_1]
  names(all_1) <- str_remove_all(names(all_1), "_1.tsv")
  all_sort <- lapply(1:length(names(all_1)), function(x) all_1[[x]][order(-all_1[[x]]$importance),])
  names(all_sort) <- names(all_1)
  return(all_sort)
}

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
      pdf(paste0(main_dir, dis_type, "/3_plots/TFs_hmp_",k, "_top_", id, ".pdf"))
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
          pdf(paste0(main_dir, dis_type, "/3_plots/TGs_hmp_",k, "_top_", id, ".pdf"),
              height=11)
        } else {
          pdf(paste0(main_dir, dis_type, "/3_plots/TGs_hmp_",k, "_top_", id, ".pdf"),
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
