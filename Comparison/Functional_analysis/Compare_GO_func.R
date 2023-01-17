# Author: Aura Zelco
# Brief description:
  # This script is used to obtain 
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all conditions
  # 3. 
  
# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # ct: celltype
  # df: dataframe


# OBS: this script is sourced in Compare_GO.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
# install.packages("tidyr")
library(clusterProfiler) # to obtain the GO information
library(org.Hs.eg.db) # the organism database to use
library(tidyr) # to re-arrange dfs
library(ggpubr) # to create a composite figure from multiple dot plots

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

# 2. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
  # Input: directory where to find ct sub-folders, if UCSC or not, list of projects ids, the individual project id to look for, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, UCSC_flag="no", individual_projs=vector(), single_proj="", ext, row_col) {
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

# 3. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, list of projects ids
  # Return: list of condition lists, each containing ct lists divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=vector()) {
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
          single_proj_list <- list(ImportCt(paste0(main_dir, folder), individual_projs = individual_projs, single_proj = single_proj))
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



# 4. Groups cts according to common annotation
# Input: list of lists generated from ImportDatasets, here combined in a vector, and the named vector used to harmonize the annotation
# Return: one dataframe containing all DEGs

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
  sex_dfs <- do.call(rbind, sex_dfs)
  return(sex_dfs)
}



# 5. Compares the GOs of a list of genes
  # Input: list of genes to be compared, which GO (BP, MO or CC)
  # Return: enriched GO (formal class compareClusterResult)

compareGO <- function(input.ls, GO.ont, gene_thresh){
  enrich <- compareCluster(geneCluster  =input.ls, 
                           fun          = "enrichGO",
                           keyType      = "SYMBOL",
                           ont          = GO.ont,
                           OrgDb        = "org.Hs.eg.db",
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  if (is.numeric(gene_thresh)) {
    enrich <- gsfilter(enrich, min=gene_thresh)
  }
  return(enrich)
}

# 6. Compares the GOs of a list of genes
  # Input: main direcotry where to save the plots, the dataframe containing all DEGs, the GO to analyze, and if a gene threshold value should be used
  # Return: nothing, saves plots and CSVs instead

CondCtSexGO <- function(main_dir, sex_df, GO.ont="BP", gene_thresh="no") {
  out_path <- paste0(main_comparison, "GO_comparison_ct_sex/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (cond in unique(sex_df$condition)) {
    print(cond)
    cond_path <- paste0(out_path, cond, "/")
    dir.create(cond_path, recursive = T, showWarnings = F)
    cond_df <- sex_df[which(sex_df$condition==cond), ]
    for (ct in unique(cond_df$common_annot)) {
      print(ct)
      f_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="F"), "gene_id"]
      m_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="M"), "gene_id"]
      sex_comp <- list("F"=f_genes, "M"=m_genes)
      try({
        print(paste0("Calculating the ", GO.ont, " results for ", ct))
        sex_ct_GO <-  compareGO(sex_comp, GO.ont, gene_thresh)
        png(paste0(cond_path, "GO_", GO.ont, "_", ct, ".png"))
        print(dotplot(sex_ct_GO, by = 'count', title=ct))
        dev.off()
        csv_ct_GO <- as.data.frame(sex_ct_GO)
        print(paste0("Saving the CSV GO results for ", ct))
        write.csv(csv_ct_GO, paste0(cond_path, "GO_", GO.ont, "_", ct, ".csv"))
      })
    }
  }
}

# 7. Creates a list, each element the gene list from each condition present for a specific ct and sex
  # Input: the dataframe, the ct and sex to be analyzed, the order in which the groups should be plotted
  # Return: the list of gene lists, from a specific ct and sex, for each condition
ExtractSexCt <- function(sex_df, ct, sex, condition_ordered) {
  sex_ct <- sex_df[which(sex_df$common_annot==ct & sex_df$sex==sex), ]
  sex_ct <- split(sex_ct$gene_id, sex_ct$condition)
  sex_order <- condition_ordered[which(condition_ordered %in% names(sex_ct))]
  sex_ct <- sex_ct[sex_order]
  return(sex_ct)
}


# 8. Compares the DEGs from all condition in a specific ct-sex group
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the GO to analyze, 
    # if a gene threshold value should be used, the order in which plot the conditions, and if the x-axis labels should be rotated by 90 degrees
  # Return: nothing, saves plots and CSVs instead

CondCtGO <- function(main_dir, sex_df, GO.ont="BP", gene_thresh="no", condition_ordered, rotate_x_axis=F) {
  out_path <- paste0(main_comparison, "GO_comparison_cts/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (ct in unique(sex_df$common_annot)) {
    print(ct)
    ct_path <- paste0(out_path, ct, "/")
    dir.create(ct_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      sex_ct <- ExtractSexCt(sex_df, ct, sex, condition_ordered)
      try({
        print(paste0("Calculating the ", GO.ont, " results for ", sex))
        sex_cond_GO <-  compareGO(sex_ct, GO.ont, gene_thresh)
        sex_cond_plot <- dotplot(sex_cond_GO, by = 'count', title=sex)
        if (rotate_x_axis) {sex_cond_plot$theme$axis.text.x$angle <- 90}
        png(paste0(ct_path, "GO_", GO.ont, "_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
        print(sex_cond_plot)
        dev.off()
        csv_cond_GO <- as.data.frame(sex_cond_GO)
        print(paste0("Saving the CSV GO results for ", sex))
        write.csv(csv_cond_GO, paste0(ct_path, "GO_", GO.ont, "_", sex, ".csv"))
      })
    }
  }
}


