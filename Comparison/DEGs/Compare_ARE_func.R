# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the AREs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
  # 1. Reads all ARE CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  
# Documentation abbreviations:
  # F and M: females and males
  # ct: celltype
  # df: dataframe
  # ds: dataset
  # ARE: androgen-response element
# OBS: this script is sourced in Compare_ARE.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(reshape) # to re-arrange the dfs
library(tidyr) # to re-arrange the dfs

# 1. Import data for each condition, asplit or not by projects
# Input: CSV files
# Return: list of dfs

ImportCondition <- function(main_dir, UCSC_flag="no", individual_projs=vector(), ext, row_col) {
  if (length(individual_projs)==0) {
    path <- paste0(main_dir, "/02B_ARE_ERE")
    str_to_remove <- "_ARE_sites.csv"
  } else {
    path <- paste0(main_dir, "/02B_ARE_ERE_proj")
    str_to_remove <- "_ARE_sites_proj.csv"
  }
  are_files <- list.files(path, full.names = F, pattern = "\\.csv$")
  are_files <- are_files[grep("ARE", are_files)]
  are_files_full <- paste(path, are_files, sep="/")
  are <- lapply(are_files_full, read.csv, row.names=1)
  names(are) <- str_remove_all(are_files, str_to_remove)
  return(are)
}

# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, list of projects ids
  # Return: list of condition lists, each containing ARE dfs divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=vector()) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      are_folder <- paste0(main_dir, folder)
    } else {
      are_folder <- paste0(main_dir, folder, "/outputs")
    }
    if (length(individual_projs)==0) {
      ds_list <- append(ds_list, list(ImportCondition(are_folder, UCSC_flag =  UCSC_flag)))
      group_names <- c(group_names, folder)
    } else {
      are_cond <- ImportCondition(are_folder, UCSC_flag =  UCSC_flag, individual_projs = individual_projs)
      are_F_split <- split(are_cond[["F"]], are_cond[["F"]]$proj)
      are_M_split <- split(are_cond[["M"]], are_cond[["M"]]$proj)
      for (proj_id in individual_projs) {
        if (proj_id %in% names(are_F_split) & proj_id %in% names(are_M_split)) {
          proj_list <- list("F"=are_F_split[[proj_id]][, -c(1,3)], "M" = are_M_split[[proj_id]][, -c(1,3)])
          ds_list <- append(ds_list, list(proj_list))
          group_names <- c(group_names, paste(folder, proj_id, sep = "_"))
        }
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}



# 3. Calculate percentages of ARE sites
  # Input: ARE filtered df
  # Return: ARE DF with updated percentages after sum due to common annotation

AREdfPerc <- function(ARE_filt) {
  ARE_filt <- transform(ARE_filt, full_perc = full * 100 / bg)
  ARE_filt <- transform(ARE_filt, half_perc = half * 100 / bg)
  ARE_filt <- transform(ARE_filt, hf_perc = hf * 100 / bg)
  ARE_filt <- transform(ARE_filt, no_overlap_perc = no_overlap * 100 / bg)
  ARE_filt_perc <- ARE_filt[, c(1, 7:10)]
  ARE_filt_perc <- melt(ARE_filt_perc, id.vars = "final")
  names(ARE_filt_perc)[names(ARE_filt_perc) == 'value'] <- 'percent'
  names(ARE_filt_perc)[names(ARE_filt_perc) == 'variable'] <- 'sites'
  levels(ARE_filt_perc$sites) <- c('Full', 'Half', 'Half-Full', 'None')
  ARE_filt_perc <- separate(ARE_filt_perc, final, into = c("ct", "condition", "sex"), remove = T, sep = "-")
  return(ARE_filt_perc)
}


# 4. Combines the datasets in one dataframe
  # Input: all datasets in one list
  # Return: all ARE dfs combined in one, averaged for common annotation

CreateAREDf <- function(ARE_all_list, common_annotation) {
  for (are_df in names(ARE_all_list)) {
    for (sex in names(ARE_all_list[[are_df]])) {
      ARE_all_list[[are_df]][[sex]]$condition <- rep(are_df, nrow(ARE_all_list[[are_df]][[sex]]))
      ARE_all_list[[are_df]][[sex]]$sex <- rep(sex, nrow(ARE_all_list[[are_df]][[sex]]))
    }
  }
  ARE_df <- do.call(rbind, unlist(ARE_all_list, recursive = F))
  rownames(ARE_df) <- NULL
  ARE_df$ct <- tolower(ARE_df$ct)
  for (ct_id in unique(ARE_df$ct)) {
    ARE_df[which(ARE_df$ct==ct_id), "ct"] <- common_annotation[ct_id] 
  }
  ARE_df$final_groups <- paste(ARE_df$ct, ARE_df$condition, ARE_df$sex, sep="-")
  ARE_filt <- data.frame()
  for (id in unique(ARE_df$final_groups)) {
    ARE_filt <- rbind(ARE_filt, c(id, colSums(ARE_df[which(ARE_df$final_groups==id), c(2:6)])))
  }
  colnames(ARE_filt) <- c("final", colnames(ARE_df)[2:6])
  ARE_filt[,2:6] <- lapply(ARE_filt[,2:6], as.numeric)
  ARE_filt <- AREdfPerc(ARE_filt)
  return(ARE_filt)
}


# 5. Plot ARE results for each ct
  # Input: ARE df for each ct
  # Return: plot

PlotARECt <- function(ARE_filt_ct, condition_ordered) {
  ARE_filt_ct$condition <- factor(ARE_filt_ct$condition, condition_ordered)
  ARE_filt_ct <- ARE_filt_ct[order(ARE_filt_ct$condition),]
  col_palette <- c("#39B600", "#9590FF","#D376FF" , "#FD61D1")
  ARE_ct_plot <- ggplot(ARE_filt_ct, aes(sex, percent, fill=sites)) +
    geom_bar(stat="identity", color="black", position = "stack") +
    #facet_wrap(~condition,scales = "free", nrow = 1) +
    facet_wrap(~condition) +
    labs(x="Developmental conditions", y="% of ARE sites", fill="Overlap ARE sites") +
    scale_fill_manual(values = c('Full' = col_palette[4] , 'Half' = col_palette[3], 'Half-Full' = col_palette[2], 'None' = col_palette[1])) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ARE_ct_plot)
}

# 6. Plot ARE results 
  # Input: ARE filtered df
  # Return: nothing, saves plot instead

PlotARE <- function(main_dir, ARE_filt, condition_ordered) {
  plot_path <- paste0(main_dir, "ARE_Across_Conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct in unique(ARE_filt$ct)) {
    print(ct)
    ct_plot <- PlotARECt(ARE_filt[which(ARE_filt$ct==ct),], condition_ordered)
    pdf(paste0(plot_path, ct, ".pdf"), width = 10)
    print(ct_plot)
    dev.off()
  } 
}


