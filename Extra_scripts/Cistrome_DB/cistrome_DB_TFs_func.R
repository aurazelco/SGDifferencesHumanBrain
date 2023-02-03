# Author: Aura Zelco
# Brief description:
  # This script is used to find the TFs overlapping with the Cistrome DB
# Brief procedure:
  # 1. Reads all DEG CSV files from the SCENIC result folder
  # 2. Intersects the TFs of each files with the reference Cistrome DB
  # 3. Saves the results as a single dataframe to the output folder
# Documentation abbreviations:
  # df: dataframe
  # db/DB: database

# OBS: this script is sourced in cistrome_DB_TFs.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(stringr) # to modify and harmonize names
library(tidyr) # to clean and re-organize dfs

# 1. Defines a custom function, used to subset using exclusion rather than inclusion criteria
`%!in%` <- Negate(`%in%`)

# 2. Import GRN files
  # Input: CSV files
  # Return: list of dfs

GRNresults <- function(in_path, res_folder) {
  if (res_folder == "1_GRN") {
    all_files <- list.files(path = in_path, pattern = "\\.tsv$",full.names = T)
    all <- lapply(all_files,function(x) {
      read.table(file = x, 
                 sep = '\t', 
                 header = TRUE)
    })
    all <- lapply(1:length(all), function(x) all[[x]][order(-all[[x]]$importance),])
    names(all) <- list.files(path = in_path, pattern = "\\.tsv$",full.names = F)
    names(all) <- str_remove_all(names(all), ".tsv")
  } else if (res_folder == "3_AUCell") {
    all_files <- list.files(path = in_path, pattern = "\\.csv$",full.names = T)
    all <- lapply(all_files,function(x) {
      read.table(file = x, 
                 sep = ',', 
                 header = TRUE)
    })
    names(all) <- list.files(path = in_path, pattern = "\\.csv$",full.names = F)
    names(all) <- str_remove_all(names(all), ".tsv")
  }
  return(all)
}

# 3. Checks which TFs from GRN output are in Cistrome DB
  # Input: list of dfs, disease flag (different file tree structure), vector of TFs to intersect our results with
  # Return: df with overlapping TFs and the names of the GRN file

CistromeFilt <- function(grn_all, dis_type=F, tf_vector) {
  grn_file <- vector()
  tfs_filt <- vector()
  for (grn_out in names(grn_all)) {
    tf_in_db <- intersect(tf_vector, grn_all[[grn_out]]$TF)
    if (length(tf_in_db) > 0) {
      grn_file <- c(grn_file, rep(grn_out, length(tf_in_db)))
      tfs_filt <- c(tfs_filt, tf_in_db)
    } else {
      grn_file <- c(grn_file, grn_out)
      tfs_filt <- c(tfs_filt, "None")
      print(paste0("The file ", grn_out, " had no TFs in the Cistrome DB and it has been assigned the value of 'None'"))
    }
  }
  cistrome_df <- data.frame(grn_file, tfs_filt)
  colnames(cistrome_df) <- c("GRN_files", "TFs_in_DB")
  if (dis_type!=F) {
    cistrome_df <- separate(cistrome_df, GRN_files, into = c("proj", "sex", "run"), sep="_", remove = F)
  } else {
    cistrome_df <- separate(cistrome_df, GRN_files, into = c("sex", "run"), sep="_", remove = F)
  }
  return(cistrome_df)
}

# 4. Check folder for Cistrome DB TFs
  # Input: list of dfs, disease flag (different file tree structure), which result folder to check, vector of TFs to intersect our results with
  # Return: Nothing - saves the dataframe to the output folder

CistromeDB <- function(main_dir, res_folder, dis_type=F, tf_vector) {
  sub_groups <- list.dirs(main_dir, full.names = F, recursive = F)
  cistrome_ls <- list()
  for (sub_group in sub_groups) {
    print(sub_group)
    in_path <- paste0(main_dir, sub_group, "/", res_folder, "/sampled_100_cells_all/")
    grn_all <- GRNresults(in_path, res_folder)
    cistrome_df <- CistromeFilt(grn_all, dis_type, tf_vector)
    out_path <- paste0(main_dir, sub_group, "/5_outputs/")
    dir.create(out_path, recursive = T, showWarnings = F)
    write.csv(cistrome_df, paste0(out_path, "Cistrome_TFs.csv"), row.names = F)
    cistrome_ls <- append(cistrome_ls, list(cistrome_df))
  }
  names(cistrome_ls) <- sub_groups
  return(cistrome_ls)
}