# Author: Aura Zelco
# Brief description:
  # This script is used to comapres the SG-DEGs with reference databases
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all groups
  # 3. Compares the DEGs with McKenzie et al 2018
  # 4. Compares the DEGs with Chlamydas et al 2022
  # 5. Compares the DEGs with SFARI

# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # ct: celltype
  # df: dataframe


# OBS: this script is sourced in Comparison_wth_ref_datasets.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

library(tidyr) # to re-arrange dfs
library(ggplot2) # to plot
library(stringr) # to format strings
library(readxl) # to read excel files
library(dplyr) # to re-arrange dfs
library(scales) # to fix palette
library(RColorBrewer) # to fix palette

# Set up the max overlap for label repelling for CNET plots
options(ggrepel.max.overlaps = Inf)


# 1. Import data for each ct
  # Input: CSV files
  # Return: list of dfs

ImportFiles <- function(path, ext, row_col) {
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
    deg <- ImportFiles(paste0(main_dir, sub_ct[ct]), ext, row_col)
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
        sex_df <- separate(sex_df, groups, into=c("groups", "sex", "ct", "gene_num"), sep="\\.")
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


# 5. Plot the presence heatmap
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotHmpRef <- function(ref_presence_df, ref_ct_id, plot_titles) {
  ref_plot <- ggplot(ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ], aes(cond_ct, gene_ids, fill=presence)) +
    geom_tile() +
    facet_grid(groups ~ sex, scales = "free") +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x="Cell types", y="Genes", fill="Genes found", title = plot_titles[ref_ct_id]) +
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
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ref_plot)
}

# 6. Plot the presence number of genes
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotBarPlotRef <- function(ref_presence_df, ref_ct_id, plot_titles) {
  ref_plot <- ggplot(ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ], aes(cond_ct, fill = presence)) +
    geom_bar(position = "stack") +
    facet_grid(sex ~ groups, scales = "free") +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x="Cell types", y="Number of Genes", fill="Genes found", title = plot_titles[ref_ct_id]) +
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
          axis.text.y = element_text(size=8, colour = "black"),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ref_plot)
}

# 7. Calculates the % of known markers in the DEGs
  # Input: the presence df, the ct to plot, the gene lists
  # Return: the percent df

RefPerc <- function(ref_presence_df, ref_ct_id, sex_df, plot_titles="no") {
  pos_markers <- ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ]
  tot_genes <- vector()
  tot_names <- vector()
  num_pos <- vector()
  for (id in unique(pos_markers$groups)) {
    for (ct in unique(pos_markers[which(pos_markers$groups==id), "cond_ct"])) {
      tot_names <- c(tot_names, paste(id, ct, "F", sep = "/"), paste(id, ct, "M", sep = "/"))
      num_pos <- c(num_pos, length(pos_markers[which(pos_markers$groups==id & pos_markers$cond_ct==ct & pos_markers$sex=="F" & pos_markers$presence=="Yes"), "gene_ids"]))
      num_pos <- c(num_pos, length(pos_markers[which(pos_markers$groups==id & pos_markers$cond_ct==ct & pos_markers$sex=="M" & pos_markers$presence=="Yes"), "gene_ids"]))
      tot_genes <- c(tot_genes, length(sex_df[which(sex_df$groups==id & sex_df$common_annot==ct & sex_df$sex=="F"), "gene_id"]))
      tot_genes <- c(tot_genes, length(sex_df[which(sex_df$groups==id & sex_df$common_annot==ct & sex_df$sex=="M"), "gene_id"]))
    }
  }
  if (plot_titles[1]=="no") {
    ref_perc <- data.frame(tot_names, num_pos, tot_genes)
  } else {
    ref_perc <- data.frame("ref_ct"=rep(plot_titles[ref_ct_id], length(tot_names)), tot_names, num_pos, tot_genes)
  }
  ref_perc <- separate(ref_perc, tot_names, into = c("groups", "ct", "sex"), sep = "/")
  ref_perc$perc <- ref_perc$num_pos * 100 / ref_perc$tot_genes
  return(ref_perc)
}

# 8. Plot the presence number of genes as percentage
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotBarPlotRefPerc <- function(ref_perc, ref_ct_id, plot_titles, groups_ordered) {
  ref_plot <- ggplot(ref_perc, aes(ct, perc, fill = factor(groups, groups_ordered[which(groups_ordered %in% groups)]))) +
    geom_bar(stat="identity", color="black") +
    facet_grid(sex ~ factor(groups, groups_ordered[which(groups_ordered %in% groups)]), scales = "free") +
    labs(x="Cell types", y="Markers %", fill="Groups", title = plot_titles[ref_ct_id]) +
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
          axis.text.y = element_text(size=8, colour = "black"),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ref_plot)
}

# 9. Plot the presence number of genes as percentage
  # Input: the presence df, the groups order
  # Return: the plot

PlotBarPlotRefPercFaceted <- function(ref_perc, groups_ordered, facets_align = "v") {
  ref_plot <- ggplot(ref_perc, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), perc, fill = ref_ct)) +
    geom_bar(stat="identity", color="black", position = position_dodge2(width = 0.9, preserve = "single")) +
    {if (facets_align=="v") facet_grid(ct ~ sex, space = "free", scales = "free")} +
    {if (facets_align=="h") facet_grid(sex ~ ct, space = "free", scales = "free")} +
    labs(x="Adult Datasets", y="Markers %", fill="Reference cell types") +
    {if (facets_align=="v") theme(panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), 
                                  panel.spacing.x=unit(0.5, "lines"),
                                  strip.text.y = element_text(size = 12, face="bold", colour = "black", angle = 0),
                                  strip.text.x = element_text(size = 12, face="bold", colour = "black"),
                                  plot.title = element_text(size=12, face="bold", colour = "black"),
                                  axis.line = element_line(colour = "black"),
                                  axis.title.x = element_text(size=12, face="bold", colour = "black"),
                                  axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
                                  axis.ticks.x=element_blank(),
                                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                                  axis.text.y = element_text(size=8, colour = "black"),
                                  axis.ticks.y = element_blank(),
                                  legend.position = "bottom", 
                                  legend.title = element_text(size=12, face="bold", colour = "black"))} +
    {if (facets_align=="h") theme(panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_blank(), 
                                  panel.spacing.x=unit(0.5, "lines"),
                                  strip.text.y = element_text(size = 12, face="bold", colour = "black", angle = 0),
                                  strip.text.x = element_text(size = 12, face="bold", colour = "black", angle = 90),
                                  plot.title = element_text(size=12, face="bold", colour = "black"),
                                  axis.line = element_line(colour = "black"),
                                  axis.title.x = element_text(size=12, face="bold", colour = "black"),
                                  axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
                                  axis.ticks.x=element_blank(),
                                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                                  axis.text.y = element_text(size=8, colour = "black"),
                                  axis.ticks.y = element_blank(),
                                  legend.position = "bottom", 
                                  legend.title = element_text(size=12, face="bold", colour = "black"))}
    
  return(ref_plot)
}


# 10. Plots if gebes from a reference df are found or not in the DEGs
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the reference df, 
    # the order in which plot the groups (and which groups to plot), the vector to use for plot titles
  # Return: nothing, saves plot instead

PlotRefCt <- function(main_dir, sex_df, ref_df, groups_ordered, ref_df_name="ref", plot_titles, facets_align){
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")
  dir.create(out_path, recursive = T, showWarnings = F)
  sex_df_filt <- sex_df[which(sex_df$groups %in% groups_ordered), ]
  presence <- vector()
  ids <- vector()
  gene_ids <- vector()
  for (sex_id in unique(sex_df_filt$sex)) {
    for (ct in unique(ref_df$Celltype)) {
      ref_genes <- ref_df[which(ref_df$Celltype==ct), "gene"]
      for (cond in unique(sex_df_filt[which(sex_df_filt$sex==sex_id), "groups"])) {
        for (ct_id in unique(sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$groups==cond), "common_annot"])) {
          presence <- c(presence, 
                        ifelse(ref_df[which(ref_df$Celltype==ct), "gene"] %in% sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$groups==cond & sex_df_filt$common_annot==ct_id), "gene_id"],
                               "Yes", "No"))
          ids <- c(ids, 
                   rep(paste(sex_id, ct, cond, ct_id, sep = "/"), length(ref_genes)))
          gene_ids <- c(gene_ids, ref_genes)
          
        }
      }
    }
  }
  ref_presence_df <- data.frame(ids, gene_ids, presence)
  ref_presence_df <- separate(ref_presence_df, ids, into=c("sex", "ref_ct", "groups", "cond_ct"), sep = "/")
  ref_presence_df$groups <- factor(ref_presence_df$groups, groups_ordered[which(groups_ordered %in% unique(ref_presence_df$groups))])
  ref_presence_df <- ref_presence_df[order(ref_presence_df$groups), ]
  #for (ref_ct_id in unique(ref_presence_df$ref_ct)) {
  #  print(plot_titles[[ref_ct_id]])
  #  pdf(paste0(out_path, plot_titles[ref_ct_id], "_hmp.pdf"), height = 15, width = 10)
  #  print(PlotHmpRef(ref_presence_df, ref_ct_id, plot_titles))
  #  dev.off()
  #  pdf(paste0(out_path, plot_titles[ref_ct_id], "_barplot.pdf"), height = 4, width = 16)
  #  print(PlotBarPlotRef(ref_presence_df, ref_ct_id, plot_titles))
  #  dev.off()
  #  ref_perc <- RefPerc(ref_presence_df, ref_ct_id, sex_df)
  #  pdf(paste0(out_path, plot_titles[ref_ct_id], "_barplot_perc.pdf"), height = 4, width = 16)
  #  print(PlotBarPlotRefPerc(ref_perc, ref_ct_id, plot_titles, groups_ordered))
  #  dev.off()
  #}
  ref_perc <- list()
  ref_perc <- lapply(1:length(unique(ref_presence_df$ref_ct)), function(x) RefPerc(ref_presence_df, unique(ref_presence_df$ref_ct)[x], sex_df, plot_titles))
  ref_perc <- do.call(rbind, ref_perc)
  ref_plot <- PlotBarPlotRefPercFaceted(ref_perc, groups_ordered, facets_align)
  if (facets_align == "v") {
    plot_param <- c(9, 14)
  } else if (facets_align=="h") {
    plot_param <- c(14, 7)
  }
  pdf(paste0(out_path, ref_df_name, "_faceted_barplot_perc_", facets_align, ".pdf"), width = plot_param[1], height = plot_param[2])
  print(ref_plot)
  dev.off()
}


# 11. Creates a dataframe with only the gens related to know diseases
  # Input: main directory where to save the file, the reference dataframe with the diseases and genes, 
    # one dataframe containing all DEGs, the reference name to be used for the output folder
  # Return: the DEG dataframe with the presence of genes-associated genes and saves the results to CSV file

CreateDisDf <- function(main_dir, ref, sex_dfs, ref_df_name) {
  dis_genes <- vector()
  group_id <- vector()
  for (dis_family in unique(ref$Disease_group)) {
    for (dis in unique(ref[which(ref$Disease_group==dis_family), "Disease"])) {
      dis_genes <- c(dis_genes, unlist(str_split(ref[which(ref$Disease_group==dis_family & ref$Disease==dis), "Affected_gene"], pattern = ", ")))
      group_id <- c(group_id, rep(paste(dis_family, dis, sep = "/"), length(unlist(str_split(ref[which(ref$Disease_group==dis_family & ref$Disease==dis), "Affected_gene"], pattern = ", ")))))
    }
  }
  ref_df <- data.frame(group_id, dis_genes)
  sex_dfs$id <- paste(sex_dfs$groups, sex_dfs$sex, sex_dfs$common_annot, sep = "/")
  dis_presence <- vector()
  dis_names <- vector()
  deg_ids <- vector()
  genes_ids <- vector()
  for (id in unique(ref_df$group_id)) {
    for (deg in unique(sex_dfs$id)) {
      genes_ids <- c(genes_ids, ref_df[which(ref_df$group_id==id), "dis_genes"])
      deg_presence <- ifelse(ref_df[which(ref_df$group_id==id), "dis_genes"] %in% sex_dfs[which(sex_dfs$id==deg), "gene_id"], "Yes", "No")
      dis_presence <- c(dis_presence, deg_presence)
      dis_names <-c(dis_names, rep(id, length(deg_presence)))
      deg_ids <- c(deg_ids, rep(deg, length(deg_presence)))
    }
  }
  ref_deg <- data.frame(dis_names, deg_ids, genes_ids, dis_presence)
  ref_deg <- separate(ref_deg, dis_names, into = c("disease_group", "disease"), sep = "/")
  ref_deg <- separate(ref_deg, deg_ids, into = c("groups", "sex", "ct"), sep = "/")
  ref_deg$dis_gene_id <- paste(ref_deg$disease, ref_deg$genes_ids, sep =  " - ")
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")  
  dir.create(out_path, showWarnings = F, recursive = T)
  write.csv(ref_deg, paste0(out_path, "disease_genes_in_degs.csv"))
  return(ref_deg)
}

# 12. Plot heatmap with the results of which disease-associated genes are found in the degs
  # Iput: the DEG data frame with the presence of genes-associated genes, the disease group to plot
  # Return: the faceted plot

PlotDisDegGroup <- function(ref_deg, dis_id, groups_ordered) {
  dis_plot <- ggplot(complete(ref_deg[which(ref_deg$disease_group==dis_id),]), aes(dis_gene_id, factor(groups, rev(groups_ordered[which(groups_ordered %in% groups)])), fill=dis_presence)) +
    geom_tile() +
    facet_grid(ct ~ sex , scales = "free") +
    labs(y="Datasets", x="Disease-associated genes", fill="Genes found", title =dis_id) +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "grey",
                      guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
}

# 13. Plots the results for the disease-associated genes for all disease groups
  # Input: main directory where to save the plots, the DEG dataframe with the presence of genes-associated genes
  # Return: nothing, saves plots instead

PlotDisDeg <- function(main_dir, ref_deg, ref_df_name, groups_ordered) {
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")  
  dir.create(out_path, showWarnings = F, recursive = T)
  for (dis in unique(ref_deg$disease_group)) {
    print(dis)
    pdf(paste0(out_path, dis, ".pdf"), width = 9, height = 14)
    print(PlotDisDegGroup(ref_deg, dis, groups_ordered))
    dev.off()
  }
}

# 14. Counts the number of repeated terms for drugs
  # Input: main directory where to save the files, the merged dataframe, the reference df, the order in which plot the groups
  # Return: count df fo SFARI genes

CountSFARI <- function(main_dir, sex_dfs, ref_df, groups_ordered) {
  path <- paste0(main_dir, "/SFARI/")
  dir.create(path, recursive = T, showWarnings = F)
  ids <- vector()
  chr_count <- vector()
  tot_degs_count <- vector()
  tot_sfari_count <- vector()
  ref_df$chr_simplified <- str_replace_all(ref_df$chr, "\\d+", "Autosome")
  for (sex in c("F", "M")) {
    for (group in unique(sex_dfs$groups)) {
      for (ct in unique(sex_dfs[which(sex_dfs$groups==group), "common_annot"])) {
        id_count <- vector()
        for (chr in unique(ref_df$chr_simplified)) {
          id_count <- c(id_count, length(intersect(ref_df[which(ref_df$chr_simplified==chr), "gene.symbol"], 
                                                   sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"])))
          ids <- c(ids, paste(group, sex, ct, chr, sep = "/"))
        }
        chr_count <- c(chr_count, id_count)
        tot_sfari_count <- c(tot_sfari_count, rep(sum(id_count), length(id_count)))
        tot_degs_count <- c(tot_degs_count, rep(length(sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"]), length(id_count)))
      }
    }
  }
  ref_count <- data.frame(ids, chr_count, tot_sfari_count, tot_degs_count)
  ref_count <- separate(ref_count, ids, into = c("groups","sex",  "ct", "chr"), sep = "/", remove = T)
  ref_count$groups <- factor(ref_count$groups, groups_ordered)
  ref_count <- ref_count[order(ref_count$groups), ]
  write.csv(ref_count, paste0(path, "SFARI_count_per_chr.csv"))
  return(ref_count)
}

# 15. Calculates hypergeometric distribution for each chr
  # Input: main directory where to save the files, count df fo SFARI genes, the background number of genes, the sfari gene counts by chromosome,
    # if split enrichment by chromosomes
  # Return: df with p-values for each combination of group/ct/sex/chr

HyperGeomSFARI <- function(main_dir, sfari_df, genes_tot, sfari_count, chr_comp = T) {
  pvalues <- vector()
  ids <- vector()
  for (group_id in unique(sfari_df$groups)) {
    for (ct_id in unique(sfari_df[which(sfari_df$groups==group_id), "ct"])) {
      for (sex_id in c("F", "M")) {
        if (chr_comp) {
          for (chr_id in names(sfari_count)) {
            pvalues <- c(pvalues, 
                         phyper(
                           sfari_df[which(sfari_df$sex==sex_id & sfari_df$chr==chr_id & sfari_df$groups==group_id & sfari_df$ct==ct_id), "chr_count"] - 1,
                           sfari_df[which(sfari_df$sex==sex_id & sfari_df$chr==chr_id & sfari_df$groups==group_id & sfari_df$ct==ct_id), "tot_degs_count"],
                           genes_tot - sfari_chr_gene_counts[[chr_id]],
                           sfari_chr_gene_counts[[chr_id]],
                           lower.tail= FALSE
                         ))
            ids <- c(ids, paste(group_id, ct_id, sex_id, chr_id, sep = "--"))
          }
        } else {
          pvalues <- c(pvalues, 
                       phyper(
                         unique(sfari_df[which(sfari_df$sex==sex_id & sfari_df$groups==group_id & sfari_df$ct==ct_id), "tot_sfari_count"]) - 1,
                         unique(sfari_df[which(sfari_df$sex==sex_id & sfari_df$groups==group_id & sfari_df$ct==ct_id), "tot_degs_count"]),
                         genes_tot - sum(sfari_count),
                         sum(sfari_count),
                         lower.tail= FALSE
                       ))
          ids <- c(ids, paste(group_id, ct_id, sex_id, sep = "--"))
        }
        
      }
    }
  }
  sfari_hypergeom <- data.frame(ids, pvalues)
  path <- paste0(main_dir, "/SFARI/")
  dir.create(path, recursive = T, showWarnings = F)
  if (chr_comp) {
    sfari_hypergeom <- separate(sfari_hypergeom, ids, into = c("groups", "ct", "sex", "chr"), sep = "--")
    write.csv(sfari_hypergeom, paste0(path, "SFARI_enrichment_pvalues_chr.csv"))
  } else {
    sfari_hypergeom <- separate(sfari_hypergeom, ids, into = c("groups", "ct", "sex"), sep = "--")
    write.csv(sfari_hypergeom, paste0(path, "SFARI_enrichment_pvalues.csv"))
  }
  sfari_hypergeom$pval_sign <- rep(NA, nrow(sfari_hypergeom))
  sfari_hypergeom[which(sfari_hypergeom$pvalues>0.05), "pval_sign"] <- "NS"
  sfari_hypergeom[which(sfari_hypergeom$pvalues<=0.05 & sfari_hypergeom$pvalues>0.01), "pval_sign"] <- "*"
  sfari_hypergeom[which(sfari_hypergeom$pvalues<=0.01 & sfari_hypergeom$pvalues>0.001), "pval_sign"] <- "**"
  sfari_hypergeom[which(sfari_hypergeom$pvalues<=0.001 & sfari_hypergeom$pvalues>0.0001), "pval_sign"] <- "***"
  sfari_hypergeom[which(sfari_hypergeom$pvalues<=0.0001), "pval_sign"] <- "****"
  sfari_hypergeom$pval_sign <- factor(sfari_hypergeom$pval_sign, c("NS","*", "**","***","****"))
  return(sfari_hypergeom)
}

# 16. Plots the SFARI count results
  # Input: main directory where to save the plots, the SFARI count df, the order in which plot the groups
  # Return: nothing, saves the plots instead

PlotSFARI <- function(main_dir, ref_count, which_comp = "tot") {
  plot_path <- paste0(main_dir, "/SFARI/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  if (which_comp=="tot") {
    pdf(paste0(plot_path, "SFARI_count_total.pdf"), width = 10, height = 15)
    print(
      ggplot(ref_count, aes(groups, tot_sfari_count, fill=sex)) +
        geom_bar(stat="identity", color="black", position = "dodge") +
        labs(x="Datasets", y="SFARI gene count", fill="") +
        facet_grid(ct ~ sex, scales = "free") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              strip.text.x = element_text(size=12, face="bold", colour = "black"),
              strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", angle = 90),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.box = "vertical",
              legend.title = element_blank())
      
    )
    dev.off()
  } else if (which_comp=="chr") {
    col_palette <- hue_pal()(3)
    pdf(paste0(plot_path, "SFARI_count_chr.pdf"), width = 10, height = 15)
    print(
      ggplot(ref_count[which(ref_count$chr_count>0), ], aes(groups, chr_count, fill=chr)) +
        geom_bar(stat="identity", color="black", position = "stack") +
        labs(x="Datasets", y="SFARI gene count", fill="Chromosomes") +
        scale_fill_manual(values = c("X" = col_palette[1],
                                     "Y"= col_palette[3],
                                     "Autosome"= col_palette[2],
                                     "X,Y"="#C77CFF")) +
        facet_grid(ct ~ sex, scales = "free") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              strip.text.x = element_text(size=12, face="bold", colour = "black"),
              strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.box = "vertical",
              legend.title = element_text(size=12, face="bold", colour = "black"))
    )
    dev.off()
  }
}

# 17. Plots the SFARI enrichment results
  # Input: main directory where to save the plots, the SFARI enriched df, 
    # the order in which plot the groups and the cell types
  # Return: nothing, saves the plots instead

PlotEnrichedPvalues <- function(main_dir, sfari_hypergeom, groups_ordered, cts_ordered, chr_comp=T) {
  plot_path <- paste0(main_dir, "/SFARI/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  brewer_palette <- brewer.pal(6,"Purples")
  if (chr_comp) {
    plot_title <- paste0(plot_path, "SFARI_hypergeom_chr.pdf")
    plt_param <- c(7,7)
  } else {
    plot_title <- paste0(plot_path, "SFARI_hypergeom.pdf")
    plt_param <- c(7,5)
  }
  pdf(plot_title, width = plt_param[1], height = plt_param[2])
  print(
    ggplot(sfari_hypergeom, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), factor(ct, rev(cts_ordered[which(cts_ordered %in% ct)])), fill=pval_sign)) +
      geom_tile(color="black") +
      {if (chr_comp) facet_grid(chr ~ sex, scales = "free")} +
      {if (chr_comp==F) facet_grid( ~ sex, scales = "free")} +
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

# 18. Counts the number of repeated terms for drugs
  # Input: main directory where to save the file, the merged dataframe, the reference df, the order in which plot the groups
  # Return: count df fo Oliva genes

CountOliva <- function(main_dir, sex_dfs, ref_df, groups_ordered, reg_split=F) {
  path <- paste0(main_dir, "/Oliva/")
  dir.create(path, recursive = T, showWarnings = F)
  ids <- vector()
  oliva_count <- vector()
  tot_degs_count <- vector()
  tot_oliva_count <- vector()
  for (sex in c("F", "M")) {
    for (group in unique(sex_dfs$groups)) {
      for (ct in unique(sex_dfs[which(sex_dfs$groups==group), "common_annot"])) {
        if (reg_split==T) {
          for (chr in unique(ref_df$chr_simplified)) {
            id_count <- vector()
            for (ctx in unique(ref_df[which(ref_df$chr_simplified==chr), "region"])) {
              ctx_chr_genes <- unique(ref_df[which(ref_df$chr_simplified==chr & ref_df$region==ctx), "HUGO_gene_id"])
              id_count <- c(id_count, length(intersect(ctx_chr_genes, 
                                                       sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"])))
              ids <- c(ids, paste(group, sex, ct, chr, ctx, sep = "/"))
            }
            oliva_count <- c(oliva_count, id_count)
            tot_oliva_count <- c(tot_oliva_count, rep(sum(id_count), length(id_count)))
            tot_degs_count <- c(tot_degs_count, rep(length(sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"]), length(id_count)))
          }
        } else {
          id_count <- vector()
          for (chr in unique(ref_df$chr_simplified)) {
            chr_genes <- unique(ref_df[which(ref_df$chr_simplified==chr), "HUGO_gene_id"])
            id_count <- c(id_count, length(intersect(chr_genes, 
                                                       sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"])))
            ids <- c(ids, paste(group, sex, ct, chr, sep = "/"))
          }
          oliva_count <- c(oliva_count, id_count)
          tot_oliva_count <- c(tot_oliva_count, rep(sum(id_count), length(id_count)))
          tot_degs_count <- c(tot_degs_count, rep(length(sex_dfs[which(sex_dfs$groups==group & sex_dfs$sex==sex & sex_dfs$common_annot==ct), "gene_id"]), length(id_count)))
        }
        
      }
    }
  }
  ref_count <- data.frame(ids, oliva_count, tot_oliva_count, tot_degs_count)
  if (reg_split==T) {
    ref_count <- separate(ref_count, ids, into = c("groups","sex",  "ct", "chr", "ctx"), sep = "/", remove = T)
    file_name <- "Oliva_count_per_chr_ctx.csv"
  } else {
    ref_count <- separate(ref_count, ids, into = c("groups","sex",  "ct", "chr"), sep = "/", remove = T)
    file_name <- "Oliva_count_per_chr.csv"
  }
  ref_count$groups <- factor(ref_count$groups, groups_ordered)
  ref_count <- ref_count[order(ref_count$groups), ]
  write.csv(ref_count, paste0(path, file_name))
  return(ref_count)
}

# 19. Calculates hypergeometric distribution for each chr
  # Input: main directory where to save the files, count df of Oliva genes, the background number of genes, the oliva gene counts by chromosome,
    # if split enrichment by chromosomes
  # Return: df with p-values for each combination of group/ct/sex/chr

HyperGeomOliva <- function(main_dir, count_df, genes_tot, ref_count, chr_comp = T) {
  pvalues <- vector()
  ids <- vector()
  for (group_id in unique(count_df$groups)) {
    for (ct_id in unique(count_df[which(count_df$groups==group_id), "ct"])) {
      for (sex_id in c("F", "M")) {
        if (chr_comp) {
          for (chr_id in names(ref_count)) {
            pvalues <- c(pvalues, 
                         phyper(
                           count_df[which(count_df$sex==sex_id & count_df$chr==chr_id & count_df$groups==group_id & count_df$ct==ct_id), "oliva_count"] - 1,
                           count_df[which(count_df$sex==sex_id & count_df$chr==chr_id & count_df$groups==group_id & count_df$ct==ct_id), "tot_degs_count"],
                           genes_tot - ref_count[[chr_id]],
                           ref_count[[chr_id]],
                           lower.tail= FALSE
                         ))
            ids <- c(ids, paste(group_id, ct_id, sex_id, chr_id, sep = "--"))
          }
        } else {
          pvalues <- c(pvalues, 
                       phyper(
                         unique(count_df[which(count_df$sex==sex_id & count_df$groups==group_id & count_df$ct==ct_id), "tot_oliva_count"]) - 1,
                         unique(count_df[which(count_df$sex==sex_id & count_df$groups==group_id & count_df$ct==ct_id), "tot_degs_count"]),
                         genes_tot - sum(ref_count),
                         sum(ref_count),
                         lower.tail= FALSE
                       ))
          ids <- c(ids, paste(group_id, ct_id, sex_id, sep = "--"))
        }
        
      }
    }
  }
  oliva_hypergeom <- data.frame(ids, pvalues)
  path <- paste0(main_dir, "Oliva/")
  dir.create(path, recursive = T, showWarnings = F)
  if (chr_comp) {
    oliva_hypergeom <- separate(oliva_hypergeom, ids, into = c("groups", "ct", "sex", "chr"), sep = "--")
    write.csv(oliva_hypergeom, paste0(path, "Oliva_enrichment_pvalues_chr.csv"))
  } else {
    oliva_hypergeom <- separate(oliva_hypergeom, ids, into = c("groups", "ct", "sex"), sep = "--")
    write.csv(oliva_hypergeom, paste0(path, "Oliva_enrichment_pvalues.csv"))
  }
  oliva_hypergeom$pval_sign <- rep(NA, nrow(oliva_hypergeom))
  oliva_hypergeom[which(oliva_hypergeom$pvalues>0.05), "pval_sign"] <- "NS"
  oliva_hypergeom[which(oliva_hypergeom$pvalues<=0.05 & oliva_hypergeom$pvalues>0.01), "pval_sign"] <- "*"
  oliva_hypergeom[which(oliva_hypergeom$pvalues<=0.01 & oliva_hypergeom$pvalues>0.001), "pval_sign"] <- "**"
  oliva_hypergeom[which(oliva_hypergeom$pvalues<=0.001 & oliva_hypergeom$pvalues>0.0001), "pval_sign"] <- "***"
  oliva_hypergeom[which(oliva_hypergeom$pvalues<=0.0001), "pval_sign"] <- "****"
  oliva_hypergeom$pval_sign <- factor(oliva_hypergeom$pval_sign, c("NS","*", "**","***","****"))
  return(oliva_hypergeom)
}

# 20. Plots the Oliva enrichment results
  # Input: main directory where to save the plots, the SFARI enriched df, 
    # the order in which plot the groups and the cell types
  # Return: nothing, saves the plots instead

PlotEnrichedPvaluesOliva <- function(main_dir, oliva_hypergeom, groups_ordered, cts_ordered, chr_comp=T) {
  plot_path <- paste0(main_dir, "Oliva/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  brewer_palette <- brewer.pal(6,"Purples")
  if (chr_comp) {
    plot_title <- paste0(plot_path, "Oliva_hypergeom_chr.pdf")
    plt_param <- c(7,7)
  } else {
    plot_title <- paste0(plot_path, "Oliva_hypergeom.pdf")
    plt_param <- c(7,5)
  }
  pdf(plot_title, width = plt_param[1], height = plt_param[2])
  print(
    ggplot(oliva_hypergeom, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), factor(ct, rev(cts_ordered[which(cts_ordered %in% ct)])), fill=pval_sign)) +
      geom_tile(color="black") +
      {if (chr_comp) facet_grid(chr ~ sex, scales = "free")} +
      {if (chr_comp==F) facet_grid( ~ sex, scales = "free")} +
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
