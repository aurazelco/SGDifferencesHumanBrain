# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease groups)
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
  # 4. Plots how many genes are found in all age groups, in all but one, etc
  # 5. Plots the total number of DEGs per ct across all groups 
  # 6. Plots the number of overlapping genes between a specific groups and all others, divided by ct and sex
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

# install.packages('heatmaply')

library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
library(reshape2) # to re-organize dfs
library(Polychrome) # to set the palette to be used in the plot


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
  # Input: directory where to find ct sub-folders, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, ext, row_col) {
  sub_ct <- list.dirs(main_dir, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste0(main_dir, sub_ct[ct]), ext, row_col)
    deg_filt <- list()
    for (k in 1:length(deg)) {
      df_filt <- as.data.frame(rownames(deg[[k]]))
      colnames(df_filt) <- c("Gene")
      deg_filt <- append(deg_filt, list(df_filt))
    }
    names(deg_filt) <- lapply(1:length(names(deg)), function(i) str_replace(names(deg)[i], "_filt", ""))
    for (i in names(deg_filt)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(deg_filt[[i]]))
        names_F <- c(names_F, sub_ct[ct])
      } else {
        df_M <- append(df_M, list(deg_filt[[i]]))
        names_M <- c(names_M, sub_ct[ct])
      }
    }
  }
  names(df_F) <- tolower(names_F)
  names(df_M) <- tolower(names_M)
  df_F <- df_F[lapply(df_F,length)>0]
  df_M <- df_M[lapply(df_M,length)>0]
  if (length(df_F) != 0 & length(df_M) != 0) {
    return(list("F"=df_F, "M"=df_M))
  } else {
    return("empty")
  }
}

# 3. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, if subfolders are present
  # Return: list of groups lists, each containing input dfs divided in F and M

ImportDatasets <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=F, ext, row_col) {
  ds_list <- list()
  ct_list <- vector()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      shared_folder <- paste0(main_dir, folder)
    } else {
      shared_folder <- paste0(main_dir, folder, "/outputs/01B_num_DEGs/")
    }
    if (individual_projs==F) {
      ds_list <- append(ds_list, list(ImportCt(shared_folder, ext, row_col)))
      ct_list <-c(ct_list, list.dirs(shared_folder, recursive=FALSE, full.names = FALSE))
      group_names <- c(group_names, folder)
    } else {
      proj_conds <- list.dirs(shared_folder, full.names = F, recursive = F)
      for (cond in proj_conds) {
        ds_list <- append(ds_list, list(ImportCt(paste0(shared_folder, "/", cond, "/outputs/01B_num_DEGs/"), ext, row_col)))
        ct_list <-c(ct_list, list.dirs(paste0(shared_folder, "/", cond, "/outputs/01B_num_DEGs/"), recursive=FALSE, full.names = FALSE))
        group_names <- c(group_names, paste(cond, folder, sep = "_"))
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(list("genes"=ds_list, "ct"=unique(ct_list)))
}


# 4. Creates the df for the input ct so that we know if a DEG is found in a certain groups or not -> used to generate hmps
  # Input: list of ct dfs, which sex and ct to analyze
  # Return: df with info whether each gene is present in all groups groups in which the ct is found

CreatePresenceCtDf <- function(sex_dfs, sex, ct) {
  sub_ct <- sex_dfs[[sex]][which(sex_dfs[[sex]]$common_annot==ct), ]
  ct_sex <-(rep(unique(sub_ct$gene_id), length(unique(sub_ct$groups))))
  ct_sex <- cbind(as.data.frame(ct_sex), rep(unique(sub_ct$groups), each=length(unique(sub_ct$gene_id))))
  colnames(ct_sex) <- c("gene_id", "groups")
  ct_sex$sex <- rep(sex, nrow(ct_sex))
  ct_sex$presence <- rep("no", nrow(ct_sex))
  for (id in unique(ct_sex$gene_id)) {
    for (groups_id in unique(ct_sex$groups)) {
      if (id %in% sub_ct[which(sub_ct$groups==groups_id), "gene_id"]) {
        ct_sex[which(ct_sex$groups==groups_id & ct_sex$gene_id==id), "presence"] <- "yes"
      }
    }
  }
  return(ct_sex)
}

# 5. Creates all PresenceDfs for all cts
  # Input: list of lists, each list corresponding to a specific groups-sex-ct combination
  # Return: list of ct dfs, with information on presence of each gene across all groups

CreatePresenceDf <- function(sex_dfs) {
  if (all(unique(sex_dfs[["F"]][,"common_annot"]) %in% unique(sex_dfs[["M"]][,"common_annot"]))) {
    ct_df_list <- list()
    for (ct in unique(sex_dfs[["F"]]$common_annot)) {
      f_ct <- CreatePresenceCtDf(sex_dfs, "F", ct)
      m_ct <- CreatePresenceCtDf(sex_dfs, "M", ct)
      df_ct <- rbind(f_ct, m_ct)
      ct_df_list <- append(ct_df_list, list(df_ct))
    }
    names(ct_df_list) <- unique(sex_dfs[["F"]]$common_annot)
    return(ct_df_list)
  } else {
    print("some cts are missing in one of the two sexes")
  }
}

# 6. Groups cts according to common annotation, then creates the presence dfs
  # Input: list of lists generated from ImportDatasets, here combined in a vector, and the named vector used to harmonize the annotation
  # Return: list of presence dfs, one per each ct

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
        sex_df <- separate(sex_df, groups, into=c("groups", "sex", "ct", "gene_num"), sep="\\.")
        sex_df$gene_num <- NULL
        sex_df$common_annot <- rep(ct, nrow(sex_df))
        sex_ct <- rbind(sex_ct, sex_df)
      }
    } 
    sex_dfs <- append(sex_dfs, list(sex_ct))
  }
  names(sex_dfs) <- c("F", "M")
  ct_df_list <- CreatePresenceDf(sex_dfs)
  return(ct_df_list)
}

# 7. Calculates the number of sex-biased DEGs in each cellular location
  # Input: list of presence dfs, one per each ct, gene location reference
  # Return: df

ExtractLocation <- function(ct_df_list, loc_ref, groups_ordered, features="no") {
  ids <- vector()
  count_df <- list()
  for (ct in names(ct_df_list)) {
    for (group in unique(ct_df_list[[ct]]$groups)) {
      for (sex in c("F", "M")) {
        if (features[1]=="no") {
          genes <- unique(ct_df_list[[ct]][which(ct_df_list[[ct]]$sex==sex & ct_df_list[[ct]]$groups==group & ct_df_list[[ct]]$presence=="yes"), "gene_id"])
        } else {
          genes <- unique(ct_df_list[[ct]][which(ct_df_list[[ct]]$sex==sex & ct_df_list[[ct]]$groups==group & ct_df_list[[ct]]$presence=="yes"), "gene_id"])
          genes <- genes[which(genes %in% features)]
        }
        sub_loc <- loc_ref[which(loc_ref$Gene %in% genes & loc_ref$Reliability!="Uncertain"), ]
        sub_loc_uncertain <- loc_ref[which(loc_ref$Gene %in% genes & loc_ref$Reliability=="Uncertain"), ]
        counts <- c(colSums(sub_loc[,4:32]), "Uncertain"=sum(colSums(sub_loc_uncertain[,4:32])))
        ids <- c(ids, paste(ct, group, sex, sep="/"))
        count_df <- append(count_df, list(counts))
      }
    }
  }
  count_df <- as.data.frame(do.call(rbind, count_df))
  count_df <- cbind(ids, count_df)
  count_df <- melt(count_df, value.name = "loc_count")
  names(count_df)[names(count_df) == 'variable'] <- 'locations'
  count_df <- separate(count_df, ids, into=c("ct", "groups", "sex"), remove = T, sep="/")
  count_df$locations <- as.character(count_df$locations)
  count_df$groups <- factor(count_df$groups, rev(groups_ordered[which(groups_ordered %in% count_df$groups)]))
  return(count_df)
}