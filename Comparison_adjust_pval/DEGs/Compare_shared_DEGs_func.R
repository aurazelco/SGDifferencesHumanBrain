# Author: Aura Zelco
  # Brief description:
    # This script is used for comparing the shared DEGs across the cts from the DEG analysis across multiple datasets (different ages/disease conditions)
  # Brief procedure:
    # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
    # 2. Manually combines the annotations to be able to compare at a general level the different groups
    # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each chromosome and sex
  # Documentation abbreviations:
    # deg: differentially expressed genes
    # F and M: females and males
    # ct: celltype
    # df: dataframe
    # ds: dataset

# OBS: this script is sourced in Compare_shared_DEGs.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs


# 1. Import data for each condition, split or not by projects
  # Input: main directory where to retrieve the files
  # Return: list of dfs

ImportCondition <- function(main_dir) {
  path <- paste0(main_dir, "/01C_num_chr")
  str_to_remove <- "_shared_genes.csv"
  shared_files <- list.files(path, full.names = F, pattern = str_to_remove)
  shared_files_full <- paste(path, shared_files, sep="/")
  shared_degs <- lapply(shared_files_full, read.csv, row.names=1)
  names(shared_degs) <- str_remove_all(shared_files, str_to_remove)
  return(shared_degs)
}

# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, if subfolders are present
  # Return: list of condition lists, each containing input dfs divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=F) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      shared_folder <- paste0(main_dir, folder)
    } else {
      shared_folder <- paste0(main_dir, folder, "/outputs")
    }
    if (individual_projs==F) {
      ds_list <- append(ds_list, list(ImportCondition(shared_folder)))
      group_names <- c(group_names, folder)
    } else {
      proj_conds <- list.dirs(shared_folder, full.names = F, recursive = F)
      for (cond in proj_conds) {
        ds_list <- append(ds_list, list(ImportCondition(paste0(shared_folder, "/", cond, "/outputs"))))
        group_names <- c(group_names, paste(cond, folder, sep = "_"))
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

# 3. Re-organizes input files into one dataframe
  # Input: list of input files
  # Return: dataframe with all the information

CreateDf <- function(shared_ls, groups_ordered) {
  all_degs <- do.call(rbind, unlist(shared_ls, recursive = F))
  ids <- as.data.frame(rownames(all_degs))
  colnames(ids) <- c("id")
  ids <- separate(ids, id, into=c("group", "sex", "num"), sep="\\.")
  all_degs <- cbind(ids[, c(1,2)], all_degs)
  rownames(all_degs) <- NULL
  all_degs$group <- factor(all_degs$group, groups_ordered)
  all_degs <- all_degs[order(all_degs$group),]
  return(all_degs)
}

# 4. Subsets the chromosome of interest and create a presence df
  # Input: dataframe with all the information, chromosome to select
  # Return: presence df

CreateChrCounts <- function(all_degs, chr_name) {
  chr <- all_degs[which(all_degs$chr_simplified==chr_name), ]
  chr_freq <- chr[!duplicated(chr[, c(1,2,4)]),c(1,2,4)]
  chr_freq <- as.data.frame(table(chr_freq$gene, dnn = c("gene")))
  chr_freq <- chr_freq[order(chr_freq$Freq, decreasing = T),]
  chr_hmp <- data.frame()
  for (id in (chr_freq$gene)) {
    sex <- unique(chr[which(chr$gene==id), "sex"])
    if (length(sex)==2) {
      sex <- "Both"
    }
    for (g in unique(chr$group)) {
      pres <- ifelse(id %in% chr[which(chr$group==g), "gene"], "Yes", "No")
      chr_hmp <- rbind(chr_hmp, 
                            c(rep(id, length(pres)),
                              rep(sex, length(pres)), 
                              rep(g, length(pres)), 
                              pres))
    }
  }
  colnames(chr_hmp) <- c("gene", "sex", "group", "presence")
  chr_hmp$gene <- factor(chr_hmp$gene, rev(chr_freq$gene))
  return(chr_hmp)
}

# 5. Plots presence df as heatmap
  # Input: main directory where to save the plots, the dataframe with all the info, the chromosome of interest, the order of the groups
  # Return: nothing, saves plot instead

PlotSharedDEGs <- function(main_dir, all_degs, chr_name, groups_ordered) {
  chr_df <- CreateChrCounts(all_degs, chr_name)
  plot_path <- paste0(main_comparison, "Shared_DEGs_across cts/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  if (chr_name=="Autosome") {
    height_param = 15
  } else {
    height_param = 5
  }
  pdf(paste0(plot_path, chr_name, ".pdf"), height = height_param)
  print(
    ggplot(chr_df, aes(factor(group, groups_ordered), gene, fill=presence)) +
      geom_tile() +
      scale_fill_manual(values = c("Yes"="#F8766D",
                                   "No"="#00BFC4"),
                        guide = guide_legend(reverse = TRUE)) +
      facet_wrap(~ factor(sex, c("F", "M", "Both")), scales = "free", nrow = 1) +
      labs(x="Groups", y="Genes", fill="Genes found") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.spacing.x=unit(0, "lines"),
            strip.text = element_text(size=12, face="bold", colour = "black"),
            plot.title = element_text(size=12, face="bold", colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=10, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.text.y = element_text(size=6, colour = "black"),
            axis.ticks.y = element_blank(),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}
