# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the AREs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
  # 1. Reads all ARE csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. Plots the percentages of ARE sites in each ct across conditions, separated by sex

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

# 1. Import data for each condition, split or not by projects
  # Input: main direcotry where to retrieve the files, if UCSC or not (different file tree structure), if ARE or ERE
  # Return: list of dfs

ImportCondition <- function(main_dir, UCSC_flag="no", ARE_ERE = ARE_ERE) {
  path <- paste0(main_dir, "/02B_ARE_ERE")
  str_to_remove <- paste0("_", ARE_ERE, "_sites.csv")
  are_ere_files <- list.files(path, full.names = F, pattern = "\\.csv$")
  are_ere_files <- are_ere_files[grep(ARE_ERE, are_ere_files)]
  are_ere_files_full <- paste(path, are_ere_files, sep="/")
  are_ere <- lapply(are_ere_files_full, read.csv, row.names=1)
  names(are_ere) <- str_remove_all(are_ere_files, str_to_remove)
  return(are_ere)
}

# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, if subfolders are present, if ARE or ERE
  # Return: list of condition lists, each containing ARE dfs divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=F, ARE_ERE) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      are_folder <- paste0(main_dir, folder)
    } else {
      are_folder <- paste0(main_dir, folder, "/outputs")
    }
    if (individual_projs==F) {
      ds_list <- append(ds_list, list(ImportCondition(are_folder, UCSC_flag =  UCSC_flag, ARE_ERE = ARE_ERE)))
      group_names <- c(group_names, folder)
    } else {
      proj_conds <- list.dirs(are_folder, full.names = F, recursive = F)
      for (cond in proj_conds) {
        ds_list <- append(ds_list, list(ImportCondition(paste0(are_folder, "/", cond, "/outputs"), UCSC_flag =  UCSC_flag, ARE_ERE = ARE_ERE)))
        group_names <- c(group_names, paste(cond, folder, sep = "_"))
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

# 3. Calculate percentages of ARE sites
  # Input: ARE filtered df
  # Return: ARE DF with updated percentages after sum due to common annotation but simplified

AREdfPercSimplified <- function(ARE_filt) {
  ARE_filt <- transform(ARE_filt, full_perc = full * 100 / bg)
  ARE_filt <- transform(ARE_filt, half_perc = half * 100 / bg)
  ARE_filt <- transform(ARE_filt, hf_perc = hf * 100 / bg)
  ARE_filt <- transform(ARE_filt, no_overlap_perc = no_overlap * 100 / bg)
  ARE_filt <- transform(ARE_filt, ARE = (full + half + hf) * 100 / bg)
  ARE_filt_perc <- ARE_filt[, c(1, 10, 11)]
  ARE_filt_perc <- melt(ARE_filt_perc, id.vars = "final")
  names(ARE_filt_perc)[names(ARE_filt_perc) == 'value'] <- 'percent'
  names(ARE_filt_perc)[names(ARE_filt_perc) == 'variable'] <- 'sites'
  levels(ARE_filt_perc$sites) <- c('None', 'ARE')
  ARE_filt_perc$sites <- factor(ARE_filt_perc$sites, c("ARE", "None"))
  ARE_filt_perc <- separate(ARE_filt_perc, final, into = c("ct", "condition", "sex"), remove = T, sep = "-")
  return(ARE_filt_perc)
}


# 4. Combines the datasets in one dataframe
  # Input: all datasets in one list
  # Return: all ARE dfs combined in one, averaged for common annotation

CreateAREDf <- function(ARE_all_list, common_annotation, simpl="no") {
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
  if (simpl=="yes") {
    ARE_filt <- AREdfPercSimplified(ARE_filt)
  } else {
    ARE_filt <- AREdfPerc(ARE_filt)
  }
  return(ARE_filt)
}

# 5. Calculate percentages of ERE sites
  # Input: ERE filtered df
  # Return: ERE DF with updated percentages after sum due to common annotation

EREdfPerc <- function(ERE_filt) {
  ERE_filt <- transform(ERE_filt, ERE_perc = ERE_overlap * 100 / bg)
  ERE_filt <- transform(ERE_filt, no_overlap_perc = no_overlap * 100 / bg)
  ERE_filt_perc <- ERE_filt[, c(1, 5:6)]
  ERE_filt_perc <- melt(ERE_filt_perc, id.vars = "final")
  names(ERE_filt_perc)[names(ERE_filt_perc) == 'value'] <- 'percent'
  names(ERE_filt_perc)[names(ERE_filt_perc) == 'variable'] <- 'sites'
  levels(ERE_filt_perc$sites) <- c('ERE', 'None')
  ERE_filt_perc <- separate(ERE_filt_perc, final, into = c("ct", "condition", "sex"), remove = T, sep = "-")
  return(ERE_filt_perc)
}

# 6. Combines the datasets in one dataframe
  # Input: all datasets in one list
  # Return: all ERE dfs combined in one, averaged for common annotation

CreateEREDf <- function(ERE_all_list, common_annotation) {
  for (ERE_df in names(ERE_all_list)) {
    for (sex in names(ERE_all_list[[ERE_df]])) {
      ERE_all_list[[ERE_df]][[sex]]$condition <- rep(ERE_df, nrow(ERE_all_list[[ERE_df]][[sex]]))
      ERE_all_list[[ERE_df]][[sex]]$sex <- rep(sex, nrow(ERE_all_list[[ERE_df]][[sex]]))
    }
  }
  ERE_df <- do.call(rbind, unlist(ERE_all_list, recursive = F))
  rownames(ERE_df) <- NULL
  ERE_df$ct <- tolower(ERE_df$ct)
  for (ct_id in unique(ERE_df$ct)) {
    ERE_df[which(ERE_df$ct==ct_id), "ct"] <- common_annotation[ct_id] 
  }
  ERE_df$final_groups <- paste(ERE_df$ct, ERE_df$condition, ERE_df$sex, sep="-")
  ERE_filt <- data.frame()
  for (id in unique(ERE_df$final_groups)) {
    ERE_filt <- rbind(ERE_filt, c(id, colSums(ERE_df[which(ERE_df$final_groups==id), c(2:4)])))
  }
  colnames(ERE_filt) <- c("final", colnames(ERE_df)[2:4])
  ERE_filt[,2:4] <- lapply(ERE_filt[,2:4], as.numeric)
  ERE_filt <- EREdfPerc(ERE_filt)
  return(ERE_filt)
}

# 7. Plot ARE results for each ct
  # Input: ARE df for each ct
  # Return: plot

PlotARECt <- function(ARE_filt_ct, groups_ordered) {
  ARE_filt_ct$condition <- factor(ARE_filt_ct$condition, groups_ordered)
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

# 8. Plot ARE results 
  # Input: main directory where to save the plots, ARE filtered df, the order of the age/condition groups
  # Return: nothing, saves plot instead

PlotARE <- function(main_dir, ARE_filt, groups_ordered) {
  plot_path <- paste0(main_dir, "ARE_ERE_Across_Conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct in unique(ARE_filt$ct)) {
    print(ct)
    ct_plot <- PlotARECt(ARE_filt[which(ARE_filt$ct==ct),], groups_ordered)
    pdf(paste0(plot_path, ct, ".pdf"))
    print(ct_plot)
    dev.off()
  } 
}

# 9. Plot ARE results faceted
  # Input: main directory where to save the plots, ARE filtered df, the order of the age/condition groups
  # Return: nothing, saves plot instead

PlotFacetedARE <- function(main_dir, ARE_filt, groups_ordered, simpl="no") {
  plot_path <- paste0(main_dir, "ARE_ERE_Across_Conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  ARE_filt$condition <- factor(ARE_filt$condition, groups_ordered)
  ARE_filt <- ARE_filt[order(ARE_filt$condition),]
  if (simpl=="no") {
    col_palette <- c("#39B600", "#9590FF","#D376FF" , "#FD61D1")
    plot_title <- "ARE_faceted.pdf"
  } else {
    col_palette <- c("#39B600", "#FD61D1")
    plot_title <- "ARE_faceted_simplified.pdf"
  }
  pdf(paste0(plot_path, plot_title), width = 20, height = 25)
  print(
    ggplot(ARE_filt, aes(sex, percent, fill=sites)) +
      geom_bar(stat="identity", color="black", position = "stack") +
      facet_grid(ct~condition, scales = "free", drop = T) +
      labs(x="Developmental conditions", y="% of ARE sites", fill="Overlap ARE sites") +
      {if (simpl=="no") scale_fill_manual(values = c('Full' = col_palette[4] , 'Half' = col_palette[3], 'Half-Full' = col_palette[2], 'None' = col_palette[1]))} +
      {if (simpl=="yes") scale_fill_manual(values = c('ARE' = col_palette[2], 'None' = col_palette[1]))} +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size=14, colour = "black", vjust = 0.7, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=14, face="bold", colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=14, face="bold", colour = "black"),
            legend.text = element_text(size=14, face="bold", colour = "black"),
            strip.text.x = element_text(size=14, face="bold", colour = "black", angle = 90),
            strip.text.y.right = element_text(size=14, face="bold", colour = "black", angle = 0))
  )
  dev.off()
}

# 10. Plot ERE results faceted
  # Input: main directory where to save the plots, ERE filtered df, the order of the age/condition groups
  # Return: nothing, saves plot instead

PlotFacetedERE <- function(main_dir, ERE_filt, groups_ordered) {
  plot_path <- paste0(main_dir, "ARE_ERE_Across_Conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  ERE_filt$condition <- factor(ERE_filt$condition, groups_ordered)
  ERE_filt <- ERE_filt[order(ERE_filt$condition),]
  col_palette <- c("#39B600", "#9590FF")  
  pdf(paste0(plot_path, "ERE_faceted.pdf"), width = 20, height = 25)
  print(
    ggplot(ERE_filt, aes(sex, percent, fill=sites)) +
      geom_bar(stat="identity", color="black", position = "stack") +
      facet_grid(ct~condition, scales = "free", drop = T) +
      labs(y="% of ERE sites", fill="Overlap ERE sites") +
      scale_fill_manual(values = c("ERE" = col_palette[2], "None" = col_palette[1])) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size=14, colour = "black", vjust = 0.7, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=14, face="bold", colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=14, face="bold", colour = "black"),
            legend.text = element_text(size=14, face="bold", colour = "black"),
            strip.text.x = element_text(size=14, face="bold", colour = "black", angle = 90),
            strip.text.y.right = element_text(size=14, face="bold", colour = "black", angle = 0))
  )
  dev.off()
}

# 11. Plot ARE and ERE results combined
  # Input: main directory where to save the plots, ARE filtered df, ERE filtered df, the order of the age/condition groups
  # Return: nothing, saves plot instead

PlotAREERECombined <- function(main_dir, ARE_filt, ERE_filt, groups_ordered, legend_cols="sites") {
  plot_path <- paste0(main_dir, "ARE_ERE_Across_Conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  ARE_filt$ARE_ERE <- rep("ARE", nrow(ARE_filt))
  ERE_filt$ARE_ERE <- rep("ERE", nrow(ERE_filt))
  comb_ARE_ERE <- rbind(ARE_filt, ERE_filt)
  if (legend_cols=="sites") {
    comb_ARE_ERE$sites <- factor(comb_ARE_ERE$sites, c("ARE", "ERE", "None"))
    pdf(paste0(plot_path, "ARE_ERE_", legend_cols ,".pdf"), width = 10, height = 12)
    print(
      ggplot(comb_ARE_ERE, aes(factor(condition, groups_ordered[which(groups_ordered %in% condition)]), percent, fill=sites)) +
        geom_bar(stat="identity", color="black", position = "stack") +
        facet_grid(ct~ARE_ERE + sex, scales = "free") +
        labs(x="Datasets", y="% of RE sites", fill="Response element sites") +
        scale_fill_manual(values = c('ARE'="#E1AD01", 'ERE'="#39B600",   "None"="#5A5A5A")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.x=element_text(size=12, face="bold", colour = "black"),
              axis.text.y = element_text(size=8, colour = "black"),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"),
              legend.text = element_text(size=12, face="bold", colour = "black"),
              strip.text.x = element_text(size=12, face="bold", colour = "black"),
              strip.text.y.right = element_text(size=12, face="bold", colour = "black", angle = 0))
    )
    dev.off()
  } else if (legend_cols=="sex"){
    comb_ARE_ERE$sex_sites <- paste(comb_ARE_ERE$sex, comb_ARE_ERE$sites, sep = " - ")
    comb_ARE_ERE$sex_sites <- factor(comb_ARE_ERE$sex_sites, 
                              c("F - ARE",  "M - ARE"  ,
                                "F - ERE",  "M - ERE",
                                "F - None", "M - None"))
    pdf(paste0(plot_path, "ARE_ERE_", legend_cols ,".pdf"), width = 10, height = 12)
    print(
      ggplot(comb_ARE_ERE, aes(factor(condition, groups_ordered[which(groups_ordered %in% condition)]), percent, fill=sex_sites)) +
        geom_bar(stat="identity", color="black", position = "stack") +
        facet_grid(ct~sex + ARE_ERE, scales = "free") +
        labs(x="Datasets", y="% of RE sites", fill="Response element sites") +
        scale_fill_manual(values = c(
          "F - ARE" ="#F8766D",  "M - ARE"="#00BFC4"  ,
          "F - ERE"= "#F8766D" ,  "M - ERE"="#00BFC4",
          "F - None" = "#5A5A5A", "M - None" = "#5A5A5A")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.x=element_text(size=12, face="bold", colour = "black"),
              axis.text.y = element_text(size=8, colour = "black"),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"),
              legend.text = element_text(size=12, face="bold", colour = "black"),
              strip.text.x = element_text(size=12, face="bold", colour = "black"),
              strip.text.y.right = element_text(size=12, face="bold", colour = "black", angle = 0))
    )
    dev.off()
  }
  
  
}