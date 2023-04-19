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
library(RColorBrewer) # to set a palette for the number of DEGs palette
library(rjson) # to import the json files containing the genes associated with the corresponding hormone
library(dplyr) # to re-organize dfs

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
# Return: list of condition lists, each containing input dfs divided in F and M

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


# 4. Groups cts according to common annotation
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

# 5. Generates a df with the nunmber of found hormone targets as absoluet numbers, percetages fo the hormone gene lists and percentages of  the sex-biased DEGs
  # Input: one dataframe with all sex-biased DEGs, the hormone gene lists in a list, and the order of the groups
  # Return: dataframe with all the counts and percentages

CreateHormonesDf <- function(sex_dfs, ref_hormones, groups_ordered) {
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
  names(hormone_ls) <- str_to_title(str_replace_all(names(ref_hormones), "/", "_"))
  hormones_df <- do.call(rbind, hormone_ls)
  hormones_df <- cbind("hormones"=gsub("\\..*", "", rownames(hormones_df)), hormones_df)
  rownames(hormones_df) <- NULL
  hormones_df <- separate(hormones_df, group_ids, into = c("condition", "ct", "sex"), sep = "/", remove = T)
  hormones_df$perc_hormones <- hormones_df$hormone_tgs * 100 / hormones_df$tot_horm
  hormones_df$perc_degs <- hormones_df$hormone_tgs * 100 / hormones_df$bg_genes
  return(hormones_df)
}

# 6. Plots the results of the hormone analysis as one faceted plot
  # Input: main directory where to save the plots, the dataframe with all the counts and percentages, 
    # the order of the groups, and the plot type between abs counts, percetages of hormones or degs
  # Return: nothing, saves the plot instead

PlotHormonesResFaceted <- function(main_dir, hormones_df, groups_ordered, plot_type="abs") {
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
    ggplot(hormones_df[which(hormones_df$yvar>0), ], aes(factor(condition, groups_ordered[which(groups_ordered %in% unique(condition))]), yvar, fill = hormones)) +
      geom_bar(stat="identity", position = "stack", color="black") +
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

# 7. Plots the results of the hormone analysis for each ct
  # Input: main directory where to save the plots, the dataframe with all the counts and percentages, 
    # the order of the groups, and the plot type between abs counts, percetages of hormones or degs
  # Return: nothing, saves the plots instead

PlotHormonesRes <- function(main_dir, hormones_df, groups_ordered, plot_type="abs") {
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
               aes(factor(condition, groups_ordered[which(groups_ordered %in% unique(condition))]), yvar, fill = ct)) +
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

# 8. Calculates the enrichment of hormones in each sub group with a hypergeometric test
  # Input: the dataframe with all the counts and percentages, the significant pvalue threshold, the minimum number of conditions to keep
  # Return: dataframe with all the significant pvalues

HormoneEnrichment <- function(hormones_df, pval_thresh=0.05, min_num_cond=1) {
  pval_names <- vector()
  pvalues <- vector()
  for (horm in unique(hormones_df$hormones)) {
    for (ct in unique(hormones_df[which(hormones_df$hormones==horm), "ct"])) {
      for (cond in unique(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct), "condition"])) {
        for (sex in c("F", "M")) {
          pval <- phyper(
            as.numeric(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct & hormones_df$condition==cond & hormones_df$sex==sex), "hormone_tgs"]) - 1,
            as.numeric(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct & hormones_df$condition==cond & hormones_df$sex==sex), "bg_genes"]),
            10000 - as.numeric(unique(hormones_df[which(hormones_df$hormones==horm), "tot_horm"])),
            as.numeric(unique(hormones_df[which(hormones_df$hormones==horm), "tot_horm"])),
            lower.tail= FALSE
          )
          pvalues <- c(pvalues, pval)
          pval_names <- c(pval_names, paste(horm, ct, cond, sex, sep="/"))
        }
      }
    }
  }
  pval_df <- data.frame(pval_names, pvalues)
  pval_df <- separate(pval_df, pval_names, into = c("hormone_id", "ct", "groups", "sex"), sep = "/", remove = T)
  pval_df<- pval_df[which(pval_df$pvalues < pval_thresh), ]
  pval_df <- pval_df %>% group_by(hormone_id) %>% filter(n() > min_num_cond)
  pval_df$pval_sign <- rep(NA, nrow( pval_df))
  pval_df[which( pval_df$pvalues>0.05), "pval_sign"] <- "NS"
  pval_df[which( pval_df$pvalues<=0.05 &  pval_df$pvalues>0.01), "pval_sign"] <- "*"
  pval_df[which( pval_df$pvalues<=0.01 &  pval_df$pvalues>0.001), "pval_sign"] <- "**"
  pval_df[which( pval_df$pvalues<=0.001 &  pval_df$pvalues>0.0001), "pval_sign"] <- "***"
  pval_df[which( pval_df$pvalues<=0.0001), "pval_sign"] <- "****"
  pval_df$pval_sign <- factor( pval_df$pval_sign, c("NS","*", "**","***","****"))
  return(as.data.frame(pval_df))
}

# 9. Plots the results of the hormone analysis as one faceted plot
  # Input: main directory where to save the plots, dataframe with all the significant pvalues, 
    # the order of the groups and the hormones to be plotted
  # Return: nothing, saves the plot instead

HmpHormoneEnrichment <- function(main_dir, pval_df, groups_ordered, features="All_hormones", plot_type="Features") {
  plot_path <- paste0(main_dir, "Hormones/Hmps/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  multi_features <- F
  if (features[1]!="All_hormones") {
    pval_df <- pval_df[which(pval_df$hormone_id %in% features),]
    plot_title <- paste0(plot_path, plot_type, ".pdf")
    if (length(features) > 1) {
      multi_features <- T
      params <- c(20, 15)
    } else {
      params <- c(10, 10)
    }
  } else {
    plot_type <- features
    plot_title <- paste0(plot_path, plot_type, ".pdf")
    params <- c(20, 15)
  }
  brewer_palette <- brewer.pal(6,"Purples")
  pdf(plot_title, height  = params[1], width = params[2])
  print(
    ggplot(complete(pval_df, hormone_id, ct, groups, sex, fill=list(pvalues=1, pval_sign="NS")), aes(factor(groups, groups_ordered[which(groups_ordered %in% unique(groups))]), ct, fill=pval_sign)) +
      geom_tile(color="black") +
      {if (features[1]=="All_hormones") facet_grid(hormone_id ~ sex)} +
      {if (features[1]!="All_hormones" & multi_features==F) facet_grid( ~ sex, scales = "free")} +
      {if (features[1]!="All_hormones" & multi_features==T) facet_grid(hormone_id ~ sex, scales = "free")} +
      scale_fill_manual(values = c("NS"="white", 
                                   "*"=brewer_palette[3],
                                   "**"=brewer_palette[4],
                                   "***"=brewer_palette[5],
                                   "****"=brewer_palette[6]),
                        na.value = "white") +
      labs(x="Datasets", y="Cell types", fill="P-values", title = str_replace_all(plot_type, "_", " ")) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            strip.text.y = element_text(size=12, colour = "black",angle = 0, face = "bold"),
            strip.text.x = element_text(size=12, colour = "black", face = "bold"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}