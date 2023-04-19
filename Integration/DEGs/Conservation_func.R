# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the results from the Conservation analysis from the DEGs workflow
# Brief procedure:
  # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Merges the result dfs in one, averaging duplicates
  # 3. Plots the resulting df
# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # ct: celltype
  # df: dataframe
  # ds: dataset

# OBS: this script is sourced in Compare_Conservation.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to re-organize the dfs


# 1. Import data for each ct
  # Input: CSV files, which conservation database, the extension
  # Return: CSV file as df

ImportConservationFile <- function(path, cons_db, threshold, ext) {
  cons_path <- paste0(path, "/", cons_db, "_fraction_in_", threshold, "_species")
  if (missing(ext)) {
    cons_file <- read.csv(paste0(cons_path, ".csv"))
  } else {
    cons_file <- read.csv(paste0(cons_path, ".", ext))
    }
  return(cons_file)
}


# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, if subfolders are present, which conservation database
  # Return: list of condition lists, each containing conservation dfs

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=F, cons_db, threshold) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      cons_folder <- paste0(main_dir, folder)
    } else {
      cons_folder <- paste0(main_dir, folder, "/outputs/02C_Conservation")
    }
    if (individual_projs==F) {
      ds_list <- append(ds_list, list(ImportConservationFile(cons_folder, cons_db, threshold)))
      group_names <- c(group_names, folder)
    } else {
      proj_conds <- list.dirs(cons_folder, full.names = F, recursive = F)
      for (cond in proj_conds) {
        ds_list <- append(ds_list, list(ImportConservationFile(paste0(cons_folder, "/", cond, "/outputs/02C_Conservation"), cons_db, threshold)))
        group_names <- c(group_names, paste(cond, folder, sep = "_"))
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

# 3. Merges all conservation dfs, retaining important metadata and averaging duplicates
  # Input: the list of conservation dfs, the common annotation to merge cell types, the order of the conditions
  # Return: averaged conservation df

CreateConservationDf <- function(cons_list, common_annotation, groups_ordered) {
  for (cond_id in names(cons_list)) {
    cons_list[[cond_id]]$condition <- rep(cond_id, nrow(cons_list[[cond_id]]))
    cons_list[[cond_id]]$common_annot <- tolower(cons_list[[cond_id]]$ct)
    for (ct_id in tolower(unique(cons_list[[cond_id]]$ct))) {
      cons_list[[cond_id]][which(tolower(cons_list[[cond_id]]$ct)==ct_id), "common_annot"] <- common_annotation[ct_id]
    }
  }
  cons_list <- do.call(rbind, cons_list)
  rownames(cons_list) <- NULL
  cons_list$grouping <- paste(cons_list$condition, cons_list$sex, cons_list$common_annot, sep="/")
  fractions_num <- vector()
  fractions_factors <- vector()
  groups <- vector()
  for (frac in c("All", "DEG_fraction")) {
    for (id in unique(cons_list$grouping)) {
      groups <- c(groups, id)
      fractions_factors <- c(fractions_factors, frac)
      fractions_num <- c(fractions_num,
                         mean(cons_list[which(cons_list$grouping==id), frac]))
    }
  }
  cons_df_mean <- data.frame("metadata"=groups, "group"=fractions_factors, "fractions"=fractions_num)
  cons_df_mean <- separate(cons_df_mean, metadata, into=c("condition", "sex", "ct"), sep="/", remove = T)
  cons_df_mean$frac_groups <- paste(cons_df_mean$sex, cons_df_mean$group, sep="_")
  cons_df_mean$condition <- factor(cons_df_mean$condition, groups_ordered)
  cons_df_mean <- cons_df_mean[order(cons_df_mean$condition), ]
  return(cons_df_mean)
}


# 4. Plots the results
  # Input: the main directory where to save the plots, averaged conservation df, which conservation database,
    # the threshold used in the DEGs analysis
  # Return: nothing, saves plot instead

PlotConservationComparison <- function(main_dir, cons_df_mean, cons_db, threshold) {
  plot_path <- paste0(main_dir, "Conservation_", cons_db, "/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  cons_df_mean <- droplevels(cons_df_mean)
  cons_df_mean$frac_groups <- str_replace_all(cons_df_mean$frac_groups, c("F_All" = "All Female Genes", "M_All" = "All Male Genes", 
                                                                          "F_DEG_fraction" = "Female-biased genes", "M_DEG_fraction" = "Male-biased genes"))
  cons_df_mean$frac_groups <- factor(cons_df_mean$frac_groups, c("All Female Genes",  "Female-biased genes", "All Male Genes", "Male-biased genes"))
  pdf(paste0(plot_path, "Conservation_comparison.pdf"), height = 15, width = 12)
  print(  
    ggplot(cons_df_mean, aes(ct, fractions * 100, fill=frac_groups)) +
            geom_bar(stat='identity', position='dodge', color="black") + 
            facet_wrap(~condition, scales = "free_x") +
            labs(x="Cell types", y="Fraction of conserved genes (%)", fill="Fraction groups") +
            #scale_fill_discrete(labels=c("All Female Genes", "Female-biased genes", "All Male Genes", "Male-biased genes")) +
            scale_fill_manual(values = c("All Female Genes" = "#D3D3D3",  
                                         "Female-biased genes" = "#F8766D", 
                                         "All Male Genes" = "#5A5A5A", 
                                         "Male-biased genes"="#00BFC4")) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  axis.title.x = element_blank(),
                  strip.text = element_text(size=12, face="bold", colour = "black"),
                  axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                  axis.ticks.x=element_blank(),
                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                  legend.position = "bottom", 
                  legend.title = element_text(size=12, face="bold", colour = "black"))
    )
  dev.off()
}