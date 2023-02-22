# Author: Aura Zelco
# Brief description:
# This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
# 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
# 2. Manually combines the annotations to be able to compare at a general level the different celltypes
# 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
# 4. Plots how many genes are found in all age groups, in all but one, etc
# 5. Plots the total number of DEGs per ct across all conditions 
# 6. Plots the number of overlapping genes between a specific condition and all others, divided by ct and sex
# Documentation abbreviations:
# deg: differentially expressed genes
# F and M: females and males
# ct: celltype
# df: dataframe
# ds: dataset
# hmps: heatmaps

# OBS: this script is sourced in Compare_DEGs.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
#library(ggpubr) # to assemble plots together before saving
#library(biomaRt) # to query to which chromosome the shared genes belong to
#library(scales) # to set the palette to be used in the PlotDEGsOverlap function
#library(RColorBrewer) # to set a palette for the number of DEGs palette
library(rjson)

# 1. Import data for each ct
  # Input: CSV files
  # Return: list of dfs

ImportDE <- function(path, ext, row_col) {
  if (missing(ext)) {
    deg_files <- list.files(path = path, pattern = "\\.csv$",full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  else {
    deg_files <- list.files(path = path, pattern = paste0("\\.",ext,"$"),full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  names_deg <- list.files(path = path, pattern = "\\.csv$",full.names = FALSE)
  names(deg) <- substr(names_deg, 1, nchar(names_deg)-4)
  return(deg)
}

# 2. Function to filter out ns genes and too low FC, and order based on FC 
  # Input: dataframe of DEGs
  # Return: gene list of significant genes as data.frame

Filter_gene <- function( order.gene.df, pval, FC) {
  logFC <- log2(1/FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  
  #If there are sig genes add index number for each sig gene
  return(data.frame("Genes"=rownames(gene.sig)))
}

# 3. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
  # Input: directory where to find ct sub-folders, if UCSC or not, list of projects ids, the individual project id to look for, 
    # the threshold for p-value and Fc if unfiltered data are imported, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, UCSC_flag="no", individual_projs=vector(), single_proj="", pval, FC, ext, row_col) {
  if (length(individual_projs)==0) {
    path <- paste0(main_dir, "/01B_num_DEGs")
  } else {
    path <- paste0(main_dir, "/01A_DEGs")
  }
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    if (length(individual_projs)>0) {
      if (length(grep(single_proj, names(deg))) > 0) {
        deg <- deg[grepl(single_proj, names(deg),fixed = T)]
      } else {
        next
      }
    }
    for (i in names(deg)) {
      if (UCSC_flag=="no") {
        if (length(individual_projs)==0) {
          colnames(deg[[i]]) <- c("Genes")
          rownames(deg[[i]]) <- NULL
          if (grepl("F", i, fixed=TRUE)){
            df_F <- append(df_F, list(deg[[i]]))
            names_F <- c(names_F, sub_ct[ct])
          } else {
            df_M <- append(df_M, list(deg[[i]]))
            names_M <- c(names_M, sub_ct[ct])
          }
        } else {
          deg_ct <- Filter_gene(deg[[i]], pval, FC)
          rownames(deg_ct) <- NULL
          if (grepl("F", i, fixed=TRUE)){
            df_F <- append(df_F, list(deg_ct))
            names_F <- c(names_F, sub_ct[ct])
          } else {
            df_M <- append(df_M, list(deg_ct))
            names_M <- c(names_M, sub_ct[ct])
          }
        }
      } else {
        deg_ct <- as.data.frame(rownames(deg[[i]]))
        colnames(deg_ct) <- c("Genes")
        rownames(deg_ct) <- NULL
        if (grepl("F", i, fixed=TRUE)){
          df_F <- append(df_F, list(deg_ct))
          names_F <- c(names_F, sub_ct[ct])
        } else {
          df_M <- append(df_M, list(deg_ct))
          names_M <- c(names_M, sub_ct[ct])
        }
      }
    }
  }
  names(df_F) <- tolower(names_F)
  names(df_M) <- tolower(names_M)
  df_F <- df_F[lapply(df_F,length)>0]
  df_M <- df_M[lapply(df_M,length)>0]
  if (length(df_F) != 0 & length(df_F) != 0) {
    return(list("F"=df_F, "M"=df_M))
  } else {
    return("empty")
  }
}

# 4. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, list of projects ids, the threshold for p-value and Fc if unfiltered data are imported
  # Return: list of condition lists, each containing ct lists divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=vector(), pval, FC) {
  ds_list <- list()
  ct_list <- vector()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      if (length(individual_projs)==0) {
        ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder))))
        ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
        group_names <- c(group_names, folder)
      } else {
        for (single_proj in individual_projs) {
          single_proj_list <- list(ImportCt(paste0(main_dir, folder), individual_projs = individual_projs, single_proj = single_proj, pval, FC))
          if (single_proj_list!="empty") {
            ds_list <- append(ds_list, single_proj_list)
            ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/01A_DEGs"), recursive=FALSE, full.names = FALSE))
            group_names <- c(group_names, paste(folder, single_proj, sep = "_"))
          }
        }
      }
    } else {
      if (length(individual_projs)==0) {
        ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder, "/outputs"), UCSC_flag)))
        ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/outputs/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
        group_names <- c(group_names, folder)
      } else {
        for (single_proj in individual_projs) {
          single_proj_list <- list(ImportCt(paste0(main_dir, folder), UCSC_flag, individual_projs = individual_projs, single_proj = single_proj))
          if (single_proj_list!="empty") {
            ds_list <- append(ds_list, single_proj_list)
            ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/outputs/01A_DEGs"), recursive=FALSE, full.names = FALSE))
            group_names <- c(group_names, paste(folder, single_proj, sep = "_"))
          }
          
        }
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(list("genes"=ds_list, "ct"=unique(ct_list)))
}

# 5. Groups cts according to common annotation
  # Input: list of lists generated from ImportDatasets, here combined in a vector, and the named vector used to harmonize the annotation
  # Return: one dataframe with all sex-biased DEGs

CreateSexDf <- function(list_ds, common_annot) {
  all <- unlist(list_ds, recursive = F)
  sex_dfs <- list()
  for (sex in c("F", "M")) {
    sex_list <- unlist(all[names(all)[which(grepl(paste0("\\.", sex), names(all)))]])
    sex_ct <- data.frame()
    for (ct in unique(common_annot)) {
      sex_filt <- list()
      for (ct_class in names(common_annot[which(common_annot==ct)])) {
        sex_filt <- append(sex_filt, sex_list[names(sex_list)[which(grepl(ct_class, names(sex_list)))]])
      }
      if (length(sex_filt)>0) {
        sex_df <- data.frame("groups" = rep(names(sex_filt), sapply(sex_filt, length)),
                             "gene_id" = unlist(sex_filt))
        rownames(sex_df) <- NULL
        sex_df <- separate(sex_df, groups, into=c("condition", "sex", "ct", "gene_num"), sep="\\.")
        sex_df$gene_num <- NULL
        sex_df$common_annot <- rep(ct, nrow(sex_df))
        sex_ct <- rbind(sex_ct, sex_df)
      }
    } 
    sex_dfs <- append(sex_dfs, list(sex_ct))
  }
  names(sex_dfs) <- c("F", "M")
  sex_dfs <- do.call(rbind, sex_dfs)
  rownames(sex_dfs) <- NULL
  return(sex_dfs)
}

CreateHormonesDf <- function(sex_dfs, ref_hormones, condition_ordered) {
  hormone_ls <- list()
  for (horm in names(ref_hormones)) {
    group_ids <- vector()
    hormone_tgs <- vector()
    bg_genes <- vector()
    for (cond in unique(sex_dfs$condition)) {
      for (ct in unique(sex_dfs[which(sex_dfs$condition==cond), "common_annot"])) {
        for (sex in c("F", "M")) {
          group_ids <- c(group_ids, paste(cond, ct, sex, sep = "/"))
          hormone_tgs <- c(hormone_tgs, 
                           length(intersect(
                             ref_hormones[[horm]], 
                             tolower(unique(sex_dfs[which(sex_dfs$condition==cond & sex_dfs$common_annot==ct &sex_dfs$sex==sex), "gene_id"])))))
          bg_genes <- c(bg_genes, length(unique(sex_dfs[which(sex_dfs$condition==cond & sex_dfs$common_annot==ct &sex_dfs$sex==sex), "gene_id"])))
        }
      }
    }
    tot_horm <- rep(length(ref_hormones[[horm]]), length(group_ids))
    horm_df <- data.frame(group_ids, hormone_tgs, tot_horm, bg_genes)
    horm_df$no_tgs <- horm_df$bg_genes - horm_df$hormone_tgs
    hormone_ls <- append(hormone_ls, list(horm_df))
  }
  names(hormone_ls) <- str_to_title(str_replace_all(names(hormones), "/", "_"))
  hormones_df <- do.call(rbind, hormone_ls)
  hormones_df <- cbind("hormones"=gsub("\\..*", "", rownames(hormones_df)), hormones_df)
  rownames(hormones_df) <- NULL
  hormones_df <- separate(hormones_df, group_ids, into = c("condition", "ct", "sex"), sep = "/", remove = T)
  hormones_df$perc_hormones <- hormones_df$hormone_tgs * 100 / hormones_df$tot_horm
  hormones_df$perc_degs <- hormones_df$hormone_tgs * 100 / hormones_df$bg_genes
  return(hormones_df)
}

PlotHormonesResFaceted <- function(main_dir, hormones_df, condition_ordered, plot_type="abs") {
  plot_path <- paste0(main_dir, "Hormones/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  if (plot_type=="abs") {
    plot_title <- "Hormones_absolute_counts.pdf"
    yvar <- "hormone_tgs"
    ytitle <- "Number of hormone genes"
  } else if (plot_type=="perc_hormones") {
    plot_title <- "Hormones_percentages.pdf"
    yvar <- "perc_hormones"
    ytitle <- "Hormone targets (%)"
  } else if (plot_type=="perc_degs") {
    plot_title <- "DEGs_percentages.pdf"
    yvar <- "perc_degs"
    ytitle <- "DEGs (%)"
  }
  hormones_df <- hormones_df[, c("hormones", "condition", "ct", "sex", yvar)]
  colnames(hormones_df) <- c("hormones", "condition", "ct", "sex", "yvar")
  pdf(paste0(plot_path, plot_title), height = 20, width = 15)
  print(
    ggplot(hormones_df[which(hormones_df$yvar>0), ], aes(factor(condition, condition_ordered[which(condition_ordered %in% unique(condition))]), yvar, fill = hormones)) +
      geom_bar(stat="identity", position = "fill", color="black") +
      facet_grid(ct ~ sex, scales = "free") +
      labs(x="Groups", y = ytitle, fill="Hormones") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}

PlotHormonesRes <- function(main_dir, hormones_df, condition_ordered, plot_type="abs") {
  if (plot_type=="abs") {
    plot_path <- paste0(main_dir, "Hormones/Absolute_counts/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    yvar <- "hormone_tgs"
    ytitle <- "Number of hormone genes"
  } else if (plot_type=="perc_hormones") {
    plot_path <- paste0(main_dir, "Hormones/Percentages_of_Hormone_targets/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    yvar <- "perc_hormones"
    ytitle <- "Hormone targets (%)"
  } else if (plot_type=="perc_degs") {
    plot_path <- paste0(main_dir, "Hormones/Percentages_of_DEGs/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    yvar <- "perc_degs"
    ytitle <- "DEGs (%)"
  }
  hormones_df <- hormones_df[, c("hormones", "condition", "ct", "sex", yvar)]
  colnames(hormones_df) <- c("hormones", "condition", "ct", "sex", "yvar")
  for (horm in unique(hormones_df$hormones)) {
    if (nrow(hormones_df[which(hormones_df$hormones==horm & hormones_df$yvar >0 ), ]) > 0 ) {
      pdf(paste0(plot_path, horm, ".pdf"), height = 20, width = 15)
      print(
        ggplot(hormones_df[which(hormones_df$hormones==horm & hormones_df$yvar > 0 ), ], 
               aes(factor(condition, condition_ordered[which(condition_ordered %in% unique(condition))]), yvar, fill = ct)) +
          geom_bar(stat="identity", position = "dodge", color="black") +
          facet_grid(ct ~ sex, scales = "free") +
          labs(x="Groups", y = ytitle, fill="Hormones", title = horm) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.title.x = element_text(size=12, face="bold", colour = "black"),
                axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    } else {
      print(paste0(horm, " had no hits in the sex-biased DEGs"))
    }
  }
}