# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the enrichment in the X and Y chromosomes from the DEGs analysis
# Brief procedure:
  # 1. Reads all enrichment csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. Plots the pvalues of the X and Y enrichments as a faceted heatmap

# Documentation abbreviations:
# F and M: females and males
# ct: celltype
# df: dataframe
# ds: dataset
# ARE: androgen-response element
# OBS: this script is sourced in Compare_XY_Enrichment.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(reshape) # to re-arrange the dfs
library(tidyr) # to re-arrange the dfs

# 1. Import data for each condition, split or not by projects
# Input: main direcotry where to retrieve the files, if UCSC or not (different file tree structure)
# Return: list of dfs

ImportCondition <- function(main_dir, UCSC_flag="no") {
  path <- paste0(main_dir, "/02A_Fisher_sex_genes/")
  if (UCSC_flag=="no") {
    str_to_remove <- "_Fisher_results.csv"
  } else {
    str_to_remove <- "_Fisher_results_v2.csv"
  }
  
  enrich_files_names <- list.files(path, full.names = F, pattern = "\\.csv$")
  enrich_files_names_full <- paste0(path, enrich_files_names)
  enrich_files <- lapply(enrich_files_names_full, read.csv, row.names=1)
  names(enrich_files) <- str_remove_all(enrich_files_names, str_to_remove)
  return(enrich_files)
}

# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
# Input: main directory, sub-folders list, if UCSC or not, if subfolders are present
# Return: list of condition lists, each containing ARE dfs divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=F) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      are_folder <- paste0(main_dir, folder)
    } else {
      are_folder <- paste0(main_dir, folder, "/outputs")
    }
    if (individual_projs==F) {
      ds_list <- append(ds_list, list(ImportCondition(are_folder, UCSC_flag =  UCSC_flag)))
      group_names <- c(group_names, folder)
    } else {
      proj_conds <- list.dirs(are_folder, full.names = F, recursive = F)
      for (cond in proj_conds) {
        ds_list <- append(ds_list, list(ImportCondition(paste0(are_folder, "/", cond, "/outputs"), UCSC_flag =  UCSC_flag)))
        group_names <- c(group_names, paste(cond, folder, sep = "_"))
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

PlotEnrichedPvalues <- function(main_dir, enrich_list, common_annot, condition_ordered) {
  enrich_df <- do.call(rbind, unlist(c(enrich_list), recursive = F, use.names = T))
  enrich_df <- cbind(str_replace_all(rownames(enrich_df), "\\.", "/"), enrich_df[, 4:5])
  rownames(enrich_df) <- NULL
  colnames(enrich_df) <- c("groups", "X", "Y")
  enrich_df <- melt(enrich_df, id.vars = "groups")
  names(enrich_df)[names(enrich_df) == 'variable'] <- 'XY'
  names(enrich_df)[names(enrich_df) == 'value'] <- 'pval'
  enrich_df <- separate(enrich_df, groups, into = c("condition", "sex", "ct"), sep = "/", remove = F)
  enrich_df$condition <- factor(enrich_df$condition, condition_ordered)
  enrich_df <- enrich_df[order(enrich_df$condition), ]
  enrich_df$common_annot <- rep(NA, nrow(enrich_df))
  for (ct_id in names(common_annot)) {
    enrich_df[which(tolower(enrich_df$ct)==ct_id), "common_annot"] <- common_annot[ct_id]
  }
  enrich_df$new_groups <- paste(enrich_df$condition, enrich_df$sex, enrich_df$common_annot, enrich_df$XY, sep = "/")
  pval_avg <- vector()
  for (id in unique(enrich_df$new_groups)) {
    pval_avg <- c(pval_avg,
                  mean(enrich_df[which(enrich_df$new_groups==id), "pval"]))
  }
  enrich_avg <- data.frame("groups"=unique(enrich_df$new_groups), "pval"=pval_avg)
  enrich_avg <- separate(enrich_avg, groups, into = c("condition", "sex", "ct", "XY"), sep = "/")
  enrich_avg$pval_bin <- ifelse(enrich_avg$pval > 0.05, "Non-significant", "Significant")
  enrich_avg$condition <- factor(enrich_avg$condition, condition_ordered)
  enrich_avg <- enrich_avg[order(enrich_avg$condition), ]
  out_path <- paste0(main_dir, "XY_enrichment/")
  dir.create(out_path, recursive = T, showWarnings = F)
  pdf(paste0(out_path, "pvalue_hmp.pdf"))
  print(
    ggplot(complete(enrich_avg, condition, sex, ct, XY), aes(condition, ct, fill=pval_bin)) +
      geom_tile() +
      facet_grid(sex ~ XY, scales = "free", drop = T) +
      scale_fill_manual(values = c("Significant"="#F8766D",
                                   "Non-significant"="#00BFC4"),
                        na.value = "grey",
                        guide = guide_legend(reverse = TRUE)) +
      labs(x="Age/Condition groups", y="Cell types", fill="P-value") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.spacing.x=unit(0, "lines"),
            plot.title = element_text(size=12, face="bold", colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
            axis.ticks.y = element_blank(),
            legend.position = "right", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}