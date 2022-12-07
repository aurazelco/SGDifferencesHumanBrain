`%!in%` <- Negate(`%in%`)

NumSexCt <- function(rds_obj) {
  if ("sex_ct" %!in% colnames(rds_obj@meta.data)) {
    rds_obj@meta.data$sex_ct <- paste(rds_obj@meta.data$sex, rds_obj@meta.data$cluster_final, sep="_")
  }
  num_sex_ct_df <- as.data.frame(table(rds_obj$sex_ct))
  num_sex_ct_df <- separate(num_sex_ct_df, Var1, into = c("sex" , "ct"), sep = "_", remove = F)
  names(num_sex_ct_df)[names(num_sex_ct_df) == 'Var1'] <- "idents"
  names(num_sex_ct_df)[names(num_sex_ct_df) == 'Freq'] <- "count"
  col_factors <- c("idents", "sex","ct")
  num_sex_ct_df[col_factors] <- lapply(num_sex_ct_df[col_factors], as.factor)  
  return(num_sex_ct_df)
}


FiltDf <- function(df, min_num_cells) {
  df <- droplevels(df)
  incomplete_ct <- vector()
  for (type in levels(df$ct)) {
    if ((nrow(subset(df, subset = ct==type))%%2!=0) | (any(subset(df, subset = ct==type)[,"count"] < min_num_cells))) {
      incomplete_ct <- c(incomplete_ct, type)
    }
  }
  df_filt <- df[df$ct %!in% incomplete_ct,]
  return(df_filt)
}

PlotFiltDf <- function(min_num, num_sex_ct_df, output_path) {
  for (min_cells in min_num) {
    num_filt <- FiltDf(num_sex_ct_df, min_cells)
    write.csv(num_filt, file = paste0(output_path, "final_filt_", min_cells, ".csv"),
              row.names = F)
    pdf(paste0(output_path, "filt_counts_", min_cells, ".pdf"), 10, 15)
    print(ggplot(num_sex_ct_df, aes(ct, count, fill=sex)) +
            geom_bar(stat="identity", position = "dodge") + 
            labs(x="", y="Nuclei count", fill="Sex") +
            geom_hline(yintercept = min_cells, linetype="dashed") +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  axis.title.x = element_text(size=12, face="bold", colour = "black"),
                  axis.text.x = element_text(size=8, colour = "black",angle = 45, vjust = 0.5, hjust=0.5),
                  axis.ticks.x=element_blank(),
                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                  legend.position = "bottom"))
    dev.off()
  }
}


CreateCellInfo <- function(rds_obj, main_dir) {
  Idents(rds_obj) <- "sex_ct"
  cell_info_rds <- data.frame()
  for (i in unique(rds_obj@meta.data$sex_ct)) {
    print(i)
    cell_id <- WhichCells(rds_obj, idents = i)
    og_group <- rep(i, length(cell_id))
    cell_info_rds <- rbind(cell_info_rds, data.frame(cell_id, og_group))
  }
  write.csv(cell_info_rds, paste0(main_dir, rds_obj@project.name, "/cell_info_", rds_obj@project.name,   ".csv"))
  return(cell_info_rds)
}

FiltTotDf <- function(df_list_dis) {
  filt_names <- vector()
  df_dis <- list()
  for (k in names(df_list_dis)) {
    if (!is.null(ncol(df_list_dis[[k]]))) {
      df_dis <- append(df_dis, list(rownames(df_list_dis[[k]][which(rowSums(as.matrix(df_list_dis[[k]]))!=0),])))
      filt_names <- c(filt_names, k)
    }
  }
  names(df_dis) <- filt_names
  cts <- vector()
  genes <- vector()
  for (id in names(df_dis)) {
    cts <- c(cts, rep(id, length(df_dis[[id]])))
    genes <- c(genes, df_dis[[id]])
  }
  tot <- data.frame(cts, genes)
  tot <- separate(tot, cts, into=c("sex", "ct"), sep ="_", remove = FALSE)
  names(tot)[names(tot) == "cts"] <- "og"
  col_factors <- c("og", "sex","ct")
  tot[col_factors] <- lapply(tot[col_factors], as.factor) 
  return(tot)
}

GenerateTotGenesDf <- function(df_list_dis, main_dir, rds_obj) {
  tot_df <- FiltTotDf(df_list_dis)
  sexes <- vector()
  cts <- vector()
  genes <- vector()
  for (sex_id in levels(tot_df$sex)) {
    for (ct_id in levels(tot_df$ct)) {
      common_genes <- tot_df[which(tot_df$sex==sex_id & tot_df$ct==ct_id), "genes"]
      genes <- c(genes, common_genes)
      sexes <- c(sexes, rep(sex_id, length(common_genes)))
      cts <- c(cts, rep(ct_id, length(common_genes)))
    }
  }
  tot_genes <- as.data.frame(cbind(sexes, cts, genes))
  colnames(tot_genes) <- c("sex", "ct", "genes")
  col_factors <- c("sex", "ct")
  tot_genes[col_factors] <- lapply(tot_genes[col_factors], as.factor) 
  write.csv(tot_genes, paste0(main_dir, rds_obj@project.name, "/tot_genes_ct_", rds_obj@project.name,   ".csv"))
}


Top2000ExprMtx <- function(mtx, main_dir, rds_obj) {
  mtx$SD <- rowSds(as.matrix(mtx))
  mtx <- mtx[which(mtx$SD > 0), ]
  if (nrow(mtx) * 0.25 > 2000) {
    mtx <- mtx[which(mtx$SD > quantile(mtx$SD)[4]), ]
    print("kept the top fourth quantile")
  } else {
    print("less than 2k genes above third quantile")
  }
  print("ordering the mtx in descending order of SD")
  mtx <- mtx[order(-mtx$SD),] 
  print("keeping only the top 2000 genes with highest SD")
  mtx <- mtx[1:2000, ]
  mtx <- cbind("Genes" = rownames(mtx), mtx)
  rownames(mtx) <- NULL
  expr_sums <- colSums(mtx[2:ncol(mtx)])
  if (identical(length(which(expr_sums>0)), length(expr_sums))) {
    print("all columns express at least one gene")
  } else {
    mtx <- mtx[ , !(names(mtx) %in% which(expr_sums>0))]
    print("calculate how many cells have been filtered out")
  }
  mtx <- mtx %>% 
    relocate(SD, .after = Genes)
  print("saving the expression mtx as RDS")
  saveRDS(mtx, paste0(main_dir, "top_2000_SD_expr_matrix_",  rds_obj@project.name, ".rds"))
  print("done")
}

RemoveDfs <- function(df_list, threshold) {
  incomplete_dfs <- vector()
  for (group_id in names(df_list)) {
    if (ncol(df_list[[group_id]]) < (threshold + 1)) {
      incomplete_dfs <- c(incomplete_dfs, group_id)
    }
  }
  add_counterpart <- vector()
  for (i in incomplete_dfs) {
    if (grepl("Female", i)) {
      m_id <- str_replace(i, "Female", "Male")
      if (m_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, m_id)
      }
    } else {
      f_id <- str_replace(i, "Male", "Female")
      if (f_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, f_id)
      }
    }
  }
  incomplete_dfs <- c(incomplete_dfs, add_counterpart)
  for (i in incomplete_dfs) {
    df_list[[i]] <- NULL
  }
  return(df_list)
}


PlotGroupNumbers <- function(df_list, thresh, main_dir) {
  df_list <- RemoveDfs(df_list, thresh)
  ids <- as.data.frame(names(df_list))
  colnames(ids) <- c("groups")
  ids <- separate(ids, groups, into = c("sex","ct"), sep="_", remove=FALSE)
  col_factors <- c("sex","ct")
  ids[col_factors] <- lapply(ids[col_factors], as.factor)  
  ids$length_groups <- sapply(1:length(names(df_list)), function(i) ncol(df_list[[i]]))
  dir.create(main_dir, showWarnings = F, recursive = T)
  pdf(paste0(main_dir, "num_filt_", thresh, "_cells.pdf"))
  print(
    ggplot(ids, aes(ct, length_groups, fill=sex)) +
      geom_bar(stat="identity", position = "dodge") + 
      geom_hline(yintercept = thresh, linetype="dashed", color = "black") +
      labs(title = paste0("Filter: ", thresh, " cells"), x="cell types", y="# of cells", fill="sex") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.5, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom",
            plot.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}

CheckDfs <- function(group_list) {
  remove_groups <- vector()
  if (length(group_list) %% 2 != 0 ) {
    f_list <- vector()
    m_list <- vector()
    for (i in names(group_list)) {
      if (grepl("Female", i)) {
        gen <- str_remove(i, "Female_")
        f_list <- c(f_list, gen)
      } else {
        gen <- str_remove(i, "Male_")
        m_list <- c(m_list, gen)
      }
    }
    if (identical(m_list, f_list) == FALSE) {
      if (length(m_list) > length(f_list)) {
        remove_groups <- c(m_list[which(m_list %!in% f_list)], "Male")
      } else {
        remove_groups <- c(f_list[which(f_list %!in% m_list)], "Female")
      }
    }
  }
  return(remove_groups)
}

RandomSampling <- function(group_list, num_sampling, num_cells, main_dir) {
  if (any(grepl("_Unknown", names(group_list)))) {
    group_list_filt <- group_list[-which(grepl("_Unknown",names(group_list)))]
  } else {
    group_list_filt <- group_list
  }
  print(names(group_list_filt))
  sampled_dfs <-list()
  sampled_names <- vector()
  for (id in names(group_list_filt)) {
    for (k in 1:num_sampling) {
      sampled <- data.frame()
      sampled <- sample(group_list_filt[[id]][-1], num_cells)
      if ("Genes" %in% colnames(group_list_filt[[id]])) {
        sampled <- cbind("Genes" = group_list_filt[[id]]$Genes, sampled)
      } else {
        sampled <- cbind("Genes" = rownames(group_list_filt[[id]]), sampled)
      }
      sampled_dfs <- append(sampled_dfs, list(sampled))
      sampled_names <- c(sampled_names, paste(id, k, sep="_"))
    }
  }
  names(sampled_dfs) <- lapply(1:length(sampled_names), function(i) str_replace_all(sampled_names[i]," ", "_"))
  dir.create(paste0(main_dir, "0_input_dfs/sampled_", num_cells, "_cells"), showWarnings = FALSE, recursive = T)
  lapply(1:length(names(sampled_dfs)), function(i) write.csv(sampled_dfs[[i]], 
                                                             file = paste0(main_dir, "0_input_dfs/sampled_", num_cells, "_cells/", names(sampled_dfs)[i], ".csv"),
                                                             row.names = FALSE))
}
