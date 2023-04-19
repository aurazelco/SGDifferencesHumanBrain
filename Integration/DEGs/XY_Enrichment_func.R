# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the enrichment in the X and Y chromosomes from the DEGs analysis
# Brief procedure:
  # 1. Reads all enrichment csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. Re-calculates enrichment for X and Y
  # 4. Plots the pvalues of the X and Y enrichments as a faceted heatmap

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
library(RColorBrewer) # for palette

# 1. Import data for each condition, split or not by projects
# Input: main direcotry where to retrieve the files, if UCSC or not (different file tree structure)
# Return: list of dfs

ImportCondition <- function(main_dir, UCSC_flag="no") {
  path <- paste0(main_dir, "/02A_HyperGeom_sex_genes/")
  if (UCSC_flag=="no") {
    str_to_remove <- "_HyperGeom_results.csv"
  } else {
    str_to_remove <- "_HyperGeom_results_v2.csv"
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

# 3. Calculate enrichment
  # Input: main directory where to save the files, list fo imported files, the background number of genes
  # Return: df with p-values for each combination of group/ct/sex/chr

HyperGeomXY <- function(main_dir, enrich_list, genes_tot, bg_num_genes, common_annot) {
  enrich_df <- do.call(rbind, unlist(c(enrich_list), recursive = F, use.names = T))
  enrich_df <- cbind(str_replace_all(rownames(enrich_df), "\\.", "/"), enrich_df[, c(1:3)])
  rownames(enrich_df) <- NULL
  colnames(enrich_df) <- c("ids", "Autosome", "X", "Y")
  enrich_df <- separate(enrich_df, ids, into=c("groups", "sex", "ct"), sep="/")
  enrich_df$common_annot <- rep(NA, nrow(enrich_df))
  for (ct_id in names(common_annot)) {
    enrich_df[which(tolower(enrich_df$ct)==ct_id), "common_annot"] <- common_annot[ct_id]
  }
  new_ids <- vector()
  auto_new <- vector()
  x_new <- vector()
  y_new <- vector()
  tot_degs <- vector()
  for (group_id in unique(enrich_df$groups)) {
    for (ct_id in unique(enrich_df[which(enrich_df$groups==group_id), "common_annot"])) {
      for (sex_id in c("F", "M")) {
        new_ids <- c(new_ids, paste(group_id, sex_id, ct_id, sep="/"))
        auto_new <- c(auto_new, sum(enrich_df[which(enrich_df$sex==sex_id & enrich_df$groups==group_id & enrich_df$common_annot==ct_id), "Autosome"]))
        x_new <- c(x_new, sum(enrich_df[which(enrich_df$sex==sex_id & enrich_df$groups==group_id & enrich_df$common_annot==ct_id), "X"]))
        y_new <- c(y_new, sum(enrich_df[which(enrich_df$sex==sex_id & enrich_df$groups==group_id & enrich_df$common_annot==ct_id), "Y"]))
      }
    }
  }
  enrich_filt <- data.frame(new_ids, auto_new, x_new, y_new)
  colnames(enrich_filt) <- c("ids", "Autosome", "X", "Y")
  enrich_filt <- separate(enrich_filt, ids, into=c("groups", "sex", "ct"), sep="/")
  enrich_filt$tot_degs <- as.integer(rowSums(enrich_filt[, c(4:6)]))
  pvalues <- vector()
  ids <- vector()
  for (group_id in unique(enrich_filt$groups)) {
    for (ct_id in unique(enrich_filt[which(enrich_filt$groups==group_id), "ct"])) {
      for (sex_id in c("F", "M")) {
        for (chr_id in names(genes_tot)) {
          pvalues <- c(pvalues, 
                       phyper(
                         enrich_filt[which(enrich_filt$sex==sex_id & enrich_filt$groups==group_id & enrich_filt$ct==ct_id), chr_id] - 1,
                         enrich_filt[which(enrich_filt$sex==sex_id & enrich_filt$groups==group_id & enrich_filt$ct==ct_id), "tot_degs"],
                         bg_num_genes - genes_tot[[chr_id]],
                         genes_tot[[chr_id]],
                         lower.tail= FALSE
                       ))
          ids <- c(ids, paste(group_id, ct_id, sex_id, chr_id, sep = "--"))
        }
      }
    }
  }
  XY_hypergeom <- data.frame(ids, pvalues)
  XY_hypergeom <- separate(XY_hypergeom, ids, into = c("groups", "ct", "sex", "chr"), sep = "--")
  path <- paste0(main_dir, "XY_enrichment//")
  dir.create(path, recursive = T, showWarnings = F)
  write.csv(XY_hypergeom, paste0(path, "XY_enrichment.csv"))
  XY_hypergeom$pval_sign <- rep(NA, nrow(XY_hypergeom))
  XY_hypergeom[which(XY_hypergeom$pvalues>0.05), "pval_sign"] <- "NS"
  XY_hypergeom[which(XY_hypergeom$pvalues<=0.05 & XY_hypergeom$pvalues>0.01), "pval_sign"] <- "*"
  XY_hypergeom[which(XY_hypergeom$pvalues<=0.01 & XY_hypergeom$pvalues>0.001), "pval_sign"] <- "**"
  XY_hypergeom[which(XY_hypergeom$pvalues<=0.001 & XY_hypergeom$pvalues>0.0001), "pval_sign"] <- "***"
  XY_hypergeom[which(XY_hypergeom$pvalues<=0.0001), "pval_sign"] <- "****"
  XY_hypergeom$pval_sign <- factor(XY_hypergeom$pval_sign, c("NS","*", "**","***","****"))
  return(XY_hypergeom)
}


# 4. Plots the p-value heatmap of XY enrichment
  # Input: main directory, the list of XY enrichment results, the common annotation and the order of the groups
  # Return: nothing, saves the plot instead

PlotEnrichedPvalues <- function(main_dir, XY_hypergeom, groups_ordered, cts_ordered) {
  out_path <- paste0(main_dir, "XY_enrichment/")
  dir.create(out_path, recursive = T, showWarnings = F)
  brewer_palette <- brewer.pal(6,"Purples")
  pdf(paste0(out_path, "pvalue_hmp.pdf"))
  print(
    ggplot(XY_hypergeom, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), factor(ct, rev(cts_ordered[which(cts_ordered %in% ct)])), fill=pval_sign)) +
      geom_tile(color="black") +
      facet_grid(chr ~ sex, scales = "free") +
      scale_fill_manual(values = c("NS"="white", 
                                   "*"=brewer_palette[3],
                                   "**"=brewer_palette[4],
                                   "***"=brewer_palette[5],
                                   "****"=brewer_palette[6]),
                        na.value = "gray") +
      labs(x="Datasets", y="Cell types", fill="P-value") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.spacing.x=unit(0.5, "lines"),
            plot.title = element_text(size=12, face="bold", colour = "black"),
            strip.text.x = element_text(size=12, face="bold", colour = "black"),
            strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
            axis.ticks.y = element_blank(),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
    
  )
  dev.off()
}