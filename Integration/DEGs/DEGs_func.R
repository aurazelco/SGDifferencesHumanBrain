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

# install.packages('heatmaply')

library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
library(ggpubr) # to assemble plots together before saving
library(biomaRt) # to query to which chromosome the shared genes belong to
library(scales) # to set the palette to be used in the PlotDEGsOverlap function
library(RColorBrewer) # to set a palette for the number of DEGs palette

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


# 4. Creates the df for the input ct so that we know if a DEG is found in a certain condition or not -> used to generate hmps
  # Input: list of ct dfs, which sex and ct to analyze
  # Return: df with info whether each gene is present in all condition groups in which the ct is found

CreatePresenceCtDf <- function(sex_dfs, sex, ct) {
  sub_ct <- sex_dfs[[sex]][which(sex_dfs[[sex]]$common_annot==ct), ]
  ct_sex <-(rep(unique(sub_ct$gene_id), length(unique(sub_ct$condition))))
  ct_sex <- cbind(as.data.frame(ct_sex), rep(unique(sub_ct$condition), each=length(unique(sub_ct$gene_id))))
  colnames(ct_sex) <- c("gene_id", "condition")
  ct_sex$sex <- rep(sex, nrow(ct_sex))
  ct_sex$presence <- rep("no", nrow(ct_sex))
  for (id in unique(ct_sex$gene_id)) {
    for (condition_id in unique(ct_sex$condition)) {
      if (id %in% sub_ct[which(sub_ct$condition==condition_id), "gene_id"]) {
        ct_sex[which(ct_sex$condition==condition_id & ct_sex$gene_id==id), "presence"] <- "yes"
      }
    }
  }
  return(ct_sex)
}

# 5. Creates all PresenceDfs for all cts
  # Input: list of lists, each list corresponding to a specific condition-sex-ct combination
  # Return: list of ct dfs, with information on presence of each gene across all conditions

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
        sex_df <- separate(sex_df, groups, into=c("condition", "sex", "ct", "gene_num"), sep="\\.")
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

# 7. Plots the heatmap for a certain ct and sex combo
  # Input: presence df, which sex to plot, and the order in which to plot the conditions
  # Return: heatmap plot

PlotDEGsConditions <- function(ct_df, ct_id, groups_ordered) {
  #ct_df <- ct_df[which(ct_df$sex==sex),]
  ct_df_ordered <- groups_ordered[which(groups_ordered %in% unique(ct_df$condition))]
  ct_df$condition <- factor(ct_df$condition, ct_df_ordered)
  ct_df <- ct_df[order(ct_df$condition), ]
  
  ct_df$condition_presence <- paste(ct_df$condition, ct_df$presence, sep="_")
  group_order <- paste(rep(sort(unique(ct_df$condition)), length(c("yes", "no"))), 
                       rep(c("yes", "no"), each=length(unique(ct_df$condition))), 
                       sep = "_")
  ct_df$condition_presence <- factor(ct_df$condition_presence, group_order)
  ct_df <- ct_df[order(ct_df$condition_presence), ] 
  gene_id_order <- unique(ct_df[which(ct_df$presence=="yes"), "gene_id"])
  ct_df$gene_id <- factor(ct_df$gene_id, gene_id_order)
  ct_df <- ct_df[order(ct_df$gene_id), ]
  ct_df$presence <- str_replace_all(ct_df$presence, c("yes"="Yes", "no"="No"))
  ct_plot <- ggplot(ct_df, aes(condition, gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x="Datasets", y="Genes", fill="Genes found", title = ct_id) +
    facet_wrap(~ sex, scales = "free") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ct_plot)
}

# 9. Generates the presence hmps for each ct, putting together F and M from the same ct and saving it as a pdf
  # Input: main directory where to save the plots, the list of presence dfs, and the order in which to plot the conditions
  # Return: nothing, saves the plot instead

PlotCts <- function(main_dir, ct_df_list, groups_ordered) {
  plot_path <- paste0(main_dir, "Hmp_DEGs_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct in names(ct_df_list)) {
    print(ct)
    #f_plot <- PlotDEGsConditions(ct_df_list[[ct]], "F", groups_ordered)
    #m_plot <- PlotDEGsConditions(ct_df_list[[ct]], "M", groups_ordered)
    #ct_plot <- ggarrange(f_plot, m_plot, common.legend = T, legend = "bottom")
    ct_plot <- PlotDEGsConditions(ct_df_list[[ct]], ct, groups_ordered)
    pdf(paste0(plot_path, ct, ".pdf"))
    print(ct_plot)
    dev.off()
  }
}

# 10. Creates the df for the input sex so that we know if a DEG is found in a certain ct or not in the condition of interest -> used to generate hmps
  # Input: list of ct dfs, which sex to analyze
  # Return: df with info whether each gene is present in the cts of the condition of interest

CreatePresenceConditionSexDf <- function(sex_dfs, sex) {
  ct_sex <-(rep(unique(sex_dfs[[sex]]$gene_id), length(unique(sex_dfs[[sex]]$common_annot))))
  ct_sex <- cbind(as.data.frame(ct_sex), rep(unique(sex_dfs[[sex]]$ct_condition), each=length(unique(sex_dfs[[sex]]$gene_id))))
  colnames(ct_sex) <- c("gene_id", "ct_condition")
  ct_sex$sex <- rep(sex, nrow(ct_sex))
  ct_sex$presence <- rep("no", nrow(ct_sex))
  for (id in unique(ct_sex$gene_id)) {
    for (condition_id in unique(ct_sex$ct_condition)) {
      if (id %in% sex_dfs[[sex]][which(sex_dfs[[sex]]$ct_condition==condition_id), "gene_id"]) {
        ct_sex[which(ct_sex$ct_condition==condition_id & ct_sex$gene_id==id), "presence"] <- "yes"
      }
    }
  }
  return(ct_sex)
}

# 10. Creates the presence dfs for both sexes, across all cts
  # Input: lists of dfs, one per sex with in each all cts
  # Return: df with the presence info

CreatePresenceConditionDf <- function(sex_dfs) {
  if (all(unique(sex_dfs[["F"]][,"common_annot"]) %in% unique(sex_dfs[["M"]][,"common_annot"]))) {
    f_cond <- CreatePresenceConditionSexDf(sex_dfs, "F")
    m_cond <- CreatePresenceConditionSexDf(sex_dfs, "M")
    return(rbind(f_cond, m_cond))
  } else {
    print("some cts are missing in one of the two sexes")
  }
}

# 11. Groups cts according to common annotation, then creates the presence dfs
  # Input: list of lists generated from ImportDatasets, here combined in a vector, and the named vector used to harmonize the annotation
  # Return: df with the presnece info for both sexes

CreateConditionDf <- function(list_ds, common_annot, condition_filt) {
  all <- unlist(list_ds[condition_filt], recursive = F)
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
        sex_df$ct_condition <- paste(sex_df$common_annot, sex_df$condition, sep="-")
        sex_ct <- rbind(sex_ct, sex_df)
      }
    } 
    sex_dfs <- append(sex_dfs, list(sex_ct))
  }
  names(sex_dfs) <- c("F", "M")
  condition_df_presence <- CreatePresenceConditionDf(sex_dfs)
  return(condition_df_presence)
}

# 12. Generate the plot for one sex of the DEGs across all cts in specific condition
  # Input: condition presence sex df
  # Return: plot

PlotAcrossConditionsSex <- function(condition_df_presence_sex, sex) {
  condition_df_presence_sex <- separate(condition_df_presence_sex, ct_condition, into = c("ct", "condition"), remove = F, sep="-")
  condition_df_presence_sex$ct_presence <- paste(condition_df_presence_sex$ct, condition_df_presence_sex$presence, sep="_")
  group_order <- paste(rep(sort(unique(condition_df_presence_sex$ct)), length(c("yes", "no"))), 
                       rep(c("yes", "no"), each=length(unique(condition_df_presence_sex$ct))), 
                       sep = "_")
  condition_df_presence_sex$ct_presence <- factor(condition_df_presence_sex$ct_presence, group_order)
  condition_df_presence_sex <- condition_df_presence_sex[order(condition_df_presence_sex$ct_presence), ] 
  gene_id_order <- unique(condition_df_presence_sex[which(condition_df_presence_sex$presence=="yes"), "gene_id"])
  condition_df_presence_sex$gene_id <- factor(condition_df_presence_sex$gene_id, gene_id_order)
  condition_df_presence_sex <- condition_df_presence_sex[order(condition_df_presence_sex$gene_id), ]
  plot_cond <- ggplot(condition_df_presence_sex, aes(ct_condition, gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("yes"="#F8766D",
                                 "no"="#00BFC4"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x="Groups", y="Genes", fill="Genes found", title = sex) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(plot_cond)
}
# 13. Generates the presence hmps for each ct, putting together F and M from the same ct and saving it as a pdf
  # Input: main directory where to save the plots, the list of presence dfs, and the order in which to plot the conditions
  # Return: nothing, saves the plot instead

PlotAcrossConditions <- function(main_dir, condition_df_presence, obj_name) {
  plot_path <- paste0(main_dir, "Hmp_DEGs_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  f_cond_plot <- PlotAcrossConditionsSex(condition_df_presence[which(condition_df_presence$sex=="F"), ], "F")
  m_cond_plot <- PlotAcrossConditionsSex(condition_df_presence[which(condition_df_presence$sex=="M"),], "M")
  condition_plot <- ggarrange(f_cond_plot, m_cond_plot, common.legend = T, legend = "bottom", nrow = 2)
  pdf(paste0(plot_path, obj_name, ".pdf"), height = 15)
  print(condition_plot)
  dev.off()
}

# 14. Creates dfs which counts in how many conditions we find each gene, per sex and ct combo
  # Input: ct df obtained previously, and the sex to analyze
  # Return: df containing for each gene the number of conditions which had that gene in their DEGs

GroupsSharingGenes <- function(ct_df, sex_id) {
  ct_df <- ct_df[which(ct_df$sex==sex_id),]
  gene_id <- vector()
  condition_count <- vector()
  for (id_gene in unique(ct_df$gene_id)) {
    gene_id <- c(gene_id, id_gene)
    condition_count <- c(condition_count, length(unique(ct_df[which(ct_df$gene_id==id_gene & ct_df$presence=="yes"), "condition"])))
  }
  sex <- rep(sex_id, length(gene_id))
  return(data.frame(gene_id, sex, condition_count))
}

# 15. Create Count Dfs for all cts
  # Input: presence df list
  # Return: list of df containing the number of conditions for each gene, for each ct

CreateCountDfs <- function(ct_df_list) {
  gene_count_dfs <- list()
  for (ct in names(ct_df_list)) {
    f_df <- GroupsSharingGenes(ct_df_list[[ct]], "F")
    m_df <- GroupsSharingGenes(ct_df_list[[ct]], "M")
    gene_count_dfs <- append(gene_count_dfs, list(rbind(f_df, m_df)))
  }
  names(gene_count_dfs) <- names(ct_df_list)
  return(gene_count_dfs)
}

# 16. Plot count dfs for each ct
# Input: shared gene df
# Return: bar plot of how many genes are shared among how many groups

PlotNumSharedGenesCt <- function(ct_df) {
  ct_plot <- ggplot(ct_df, aes(condition_count, fill=as.factor(condition_count))) +
    geom_bar(color="black") +
    facet_wrap(~ sex, scales = "free") +
    scale_x_continuous(breaks=seq(min(ct_df$condition_count), max(ct_df$condition_count),by=1)) +
    labs(y="Gene absolute count", fill="Number of datasets sharing genes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ct_plot)
}

# 17. Plot count dfs for each ct
  # Input: shared gene df
  # Return: bar plot of how many genes are shared among how many groups

PlotNumSharedGenesTot <- function(shared_genes_chr) {
  tot_plot <- ggplot(shared_genes_chr, aes(condition_count, fill=as.factor(condition_count))) +
    geom_bar(color="black") +
    facet_grid(ct ~ sex, scales = "free") +
    scale_x_continuous(breaks=seq(min(shared_genes_chr$condition_count), max(shared_genes_chr$condition_count),by=1)) +
    labs(y="Gene absolute count", fill="Number of conditions sharing genes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text = element_text(size=8, colour = "black", face="bold"),
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(tot_plot)
}

# 18. Plot fractions of shared chromosomes
  # Input: shared gene df, the color palette to use
  # Return: bar plot with fractions of chr genes shared in how many groups

PlotNumSharedGenesChr <- function(shared_genes_chr, col_palette) {
  chr_plot <- ggplot(shared_genes_chr, aes(condition_count, fill=chr_simplified)) +
    geom_bar(position = "stack", color="black") +
    scale_y_log10() +
    facet_grid(ct ~ sex, scales = "free") +
    scale_fill_manual(values = c("X" = col_palette[1],
                                 "Y"= col_palette[3],
                                 "Autosome"= col_palette[2])) +
    scale_x_continuous(breaks=seq(min(shared_genes_chr$condition_count), max(shared_genes_chr$condition_count), by=1)) +
    labs(x = "Number of datasets sharing genes", y="Log10 counts of shared genes", fill="Chromosomes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0.5, "lines"),
          strip.text.x = element_text(size=12, colour = "black", face="bold"),
          strip.text.y = element_text(size=12, colour = "black", face="bold", angle=0),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(chr_plot)
}


# 19. Plot all count dfs for all cts
  # Input: main directory where to save the plots, the list of count dfs
  # Return: nothing, saves the plot instead

PlotNumSharedGenes <- function(main_dir, gene_count_dfs) {
  plot_path <- paste0(main_dir, "Num_Shared_DEGs_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  chr_results <- lapply(gene_count_dfs, function(x) Annot.chr.name(x$gene_id))
  shared_genes_chr <- list()
  for(i in names(gene_count_dfs)){
    shared_genes_ct <- map_chr(gene_count_dfs[[i]], chr_results[[i]])
    shared_genes_chr <- append(shared_genes_chr, list(shared_genes_ct))
  }
  names(shared_genes_chr) <- names(gene_count_dfs)
  shared_genes_chr <- do.call(rbind, shared_genes_chr)
  shared_genes_chr$chr_simplified <- shared_genes_chr$chromosome_name
  shared_genes_chr$chr_simplified[which(shared_genes_chr$chr_simplified!= "X" & shared_genes_chr$chr_simplified!= "Y")] <- "Autosome"
  shared_genes_chr <- cbind("ct" = gsub('\\..*', '', rownames(shared_genes_chr)), shared_genes_chr)
  rownames(shared_genes_chr) <- NULL
  col_palette <- hue_pal()(3)
  pdf(paste0(plot_path, "tot_shared_genes.pdf"), height = 15, width = 10)
  print(PlotNumSharedGenesTot(shared_genes_chr))
  dev.off()
  pdf(paste0(plot_path, "chr_shared_genes.pdf"), height = 15, width = 10)
  print(PlotNumSharedGenesChr(shared_genes_chr, col_palette))
  dev.off()
}

# 20. Extract the genes shared among a minimum percentage of conditions
  # Input: the list of count dfs
  # Return: the list of filtered genes (of which cts that are shared across more than 1 condition)

ExtractSharedGenes <- function(gene_count_dfs, min_sharing=0.75, min_num_cond=1) {
  gene_count_filt <- list()
  for (ct in names(gene_count_dfs)) {
    if (max(gene_count_dfs[[ct]]$condition_count) > min_num_cond) {
      shared_genes <- gene_count_dfs[[ct]][which(gene_count_dfs[[ct]]$condition_count >= ceiling(max(gene_count_dfs[[ct]]$condition_count) * min_sharing)), c( "sex", "gene_id")]
      gene_count_filt <- append(gene_count_filt, list(cbind("ct"=rep(ct, nrow(shared_genes)), shared_genes)))
    }
  }
  return(gene_count_filt)
}


# 21. Function to get chromosome number from gene symbol
  # Input: the genes as vector
  # Return: the annotated genes

Annot.chr.name <- function(gene.list){
  # define biomart object
  mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "www")
  Annot_idf <- getBM(attributes = c("hgnc_symbol",
                                    "chromosome_name"),
                     filters = c("hgnc_symbol") ,
                     values=list(gene.list),
                     mart = mart)
  #delete chromosome name with CHR label
  Annot_df <- Annot_idf[!str_detect(Annot_idf$chromosome_name,  "CHR"),]
  return(Annot_df)
}

# 22. Map genes from intersected genes against chromosome
  # Input: the list of filtered genes, the annotated genes
  # Return: merged dataframe

map_chr <- function(gene_count_filt, Annot_df){
  map_chr_df <- merge(gene_count_filt, Annot_df, by.x= "gene_id", by.y= "hgnc_symbol")
  return(map_chr_df)
}

# 23. Saves the annotated shared genes to a CSV file 
  # Input: main directory where to save the CSV file, the list of filtered genes
  # Return: nothing, saves the CSV file instead

SaveSharedGenes <- function(main_dir, gene_count_dfs, min_sharing=0.75, min_num_cond=1) {
  out_path <- paste0(main_dir, "Num_Shared_DEGs_across_conditions/")
  dir.create(out_path, showWarnings = F, recursive = T)
  shared_genes <- ExtractSharedGenes(gene_count_dfs, min_sharing, min_num_cond)
  chr_results <- lapply(shared_genes, function(x) Annot.chr.name(x$gene_id))
  shared_genes_chr <- list()
  for(i in 1:length(shared_genes)){
    shared_genes_ct <- map_chr(shared_genes[[i]], chr_results[[i]])
    shared_genes_chr <- append(shared_genes_chr, list(shared_genes_ct))
  }
  shared_genes_chr <- do.call(rbind, shared_genes_chr)
  shared_genes_chr$chr_simplified <- shared_genes_chr$chromosome_name
  shared_genes_chr$chr_simplified[which(shared_genes_chr$chr_simplified!= "X" & shared_genes_chr$chr_simplified!= "Y")] <- "Autosome"
  shared_genes_chr$num_conditions <- rep(NA, nrow(shared_genes_chr))
  for (ct_id in names(gene_count_dfs)) {
      shared_genes_chr[which(shared_genes_chr$ct==ct_id), "num_conditions"] <- max(gene_count_dfs[[ct_id]]$condition_count)
  }
  write.csv(shared_genes_chr, paste0(out_path, "Shared_genes.csv"))
  #return(shared_genes_chr)
}


# 24. Count DEGs for each ct in each age
  # Input: list of presence dfs, one per each ct, order of the condition
  # Return: dataframe with num of DEGs for each ct and condition

NumDEGsAcrossConditions <- function(ct_df_list, groups_ordered) {
  ct <- vector()
  condition <- vector()
  count_degs <- vector()
  for (ct_id in names(ct_df_list)) {
    sub_ct <- ct_df_list[[ct_id]]
    for (cond in groups_ordered) {
      ct <- c(ct, rep(ct_id, 2))
      condition <- c(condition, rep(cond, 2))
      if (cond %in% unique(sub_ct$condition)) {
        count_degs <- c(count_degs, nrow(sub_ct[which(sub_ct$condition==cond & sub_ct$sex=="F" & sub_ct$presence=="yes"), ]))
        count_degs <- c(count_degs, nrow(sub_ct[which(sub_ct$condition==cond & sub_ct$sex=="M" & sub_ct$presence=="yes"), ]))
      } else {
        count_degs <- c(count_degs, rep(0, 2))
      }
    }
  }
  sex <- rep(c("F", "M"), length(condition)/2)
  num_deg_df <- data.frame(ct, condition, sex, count_degs)
  groups_ordered <- groups_ordered[which(groups_ordered %in% unique(num_deg_df$condition))]
  num_deg_df$condition <- factor(num_deg_df$condition, groups_ordered)
  num_deg_df <- num_deg_df[order(num_deg_df$condition), ]
  return(num_deg_df)
}


# 25. Plot the total num of DEGs per ct
  # Input: dataframe with num of DEGs of one ct
  # Return: plot

PlotNumDEGsCt <- function(ct_degs) {
  ct_deg_plot <- ggplot(ct_degs, aes(condition, count_degs, fill=condition)) +
    geom_bar(stat = "identity") +
    labs(y="Number of DEGs", fill="Conditions") +
    facet_wrap(~sex, scales = "free") +
    theme(
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = NA, color = "black"), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=14, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=14, face="bold", colour = "black"),
          axis.text.y = element_text(size=10, colour = "black", vjust = 0.7, hjust=0.5),
          axis.title.x = element_text(size=14, face="bold", colour = "black"),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text = element_text(size = 10),
          legend.position = "bottom", 
          legend.text = element_text(size=10, colour = "black"),
          legend.title = element_text(size=14, face="bold", colour = "black"))
  return(ct_deg_plot)
}

# 26. Plot the total num of DEGs per ct for all cts
  # Input: main directory where to save the plots, the dataframe with num of DEGs for each ct and condition
  # Return: nothing, saves the plot instead

PlotNumDEGs <- function(main_dir, num_degs_ct) {
  plot_path <- paste0(main_dir, "Num_Total_DEGs_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct_id in unique(num_degs_ct$ct)) {
    print(ct_id)
    pdf(paste0(plot_path, ct_id, ".pdf"))
    print(PlotNumDEGsCt(num_degs_ct[which(num_degs_ct$ct==ct_id),]))
    dev.off()
  }
}

# 27. Plots the total number of DEGs across conditions, faceted by ct and sex
  # Input: main directory where to save the plots, the dataframe with num of DEGs for each ct and condition
  # Return: nothing, saves the plot instead

PlotNumDEGsFacets <- function(main_dir, num_degs_ct, col_palette) {
  plot_path <- paste0(main_dir, "Num_Total_DEGs_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  pdf(paste0(plot_path, "Faceted_tot_DEGs.pdf"), height = 22, width = 15)
  print(
    ggplot(num_degs_ct, aes(condition, count_degs, fill=condition)) +
      geom_bar(stat = "identity", show.legend = T, color="black") +
      labs(x="", y="Number of DEGs", fill="Datasets") +
      scale_fill_manual(values = col_palette) + 
      facet_grid(ct ~ sex, scales = "free", switch = "y", drop = T) +
      theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = NA, color = "black"), 
            panel.spacing.x = unit(0.5, "lines"),
            strip.text = element_text(size=14, face="bold", colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.title.y = element_text(size=14, face="bold", colour = "black"),
            axis.text.y = element_text(size=10, colour = "black", vjust = 0.7, hjust=0.5),
            axis.title.x = element_text(size=14, face="bold", colour = "black"),
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.position = "bottom", 
            legend.text = element_text(size=10, colour = "black"),
            legend.title = element_text(size=14, face="bold", colour = "black"))
  )
  dev.off()
}

# 28. Plot the num of shared genes between one condition and all others
  # Input: main directory where to save the plots, list of presence dfs, one per each ct, and the order in which to plot the conditions
  # Return: nothing, saves plot and CSV instead

PlotDEGsOverlap <- function(main_dir, ct_df_list, groups_ordered) {
  plot_path <- paste0(main_dir, "DEGs_Overlap_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  filt_df <- do.call(rbind, ct_df_list)
  filt_df$presence <- ifelse(filt_df$presence=="yes", 1, 0)
  filt_df$ct <- gsub('\\..*', '', rownames(filt_df))
  for (sex_id in c("F", "M")) {
    print(sex_id)
    sex_df <- subset(filt_df, sex==sex_id & presence==1)
    ct_vec <- vector()
    comp_list <- list()
    ct_names <- vector()
    for (ct in unique(sex_df$ct)) {
      if (length(unique(sex_df[which(sex_df$ct==ct), "condition"]))>1) {
        ct_names <- c(ct_names, ct)
        print(ct)
        ct_cond <- groups_ordered[which(groups_ordered %in% unique(sex_df[which(sex_df$ct==ct), "condition"]))]
        ct_df <- data.frame()
        for (cond in ct_cond) {
          ref_genes <- sex_df[which(sex_df$ct==ct & sex_df$condition==cond), "gene_id"]
          comp_vec <- vector()
          num_common_genes <- vector()
          for (other_cond  in ct_cond[!ct_cond == cond]) {
            comp_vec <- c(comp_vec, paste(cond, other_cond, sep = " - "))
            num_common_genes <- c(num_common_genes, length(intersect(ref_genes, sex_df[which(sex_df$ct==ct & sex_df$condition==other_cond), "gene_id"])))
          }
          ct_df <- rbind(ct_df, data.frame(comp_vec, num_common_genes))
        }
        comp_list <- append(comp_list, list(ct_df))
      }
    }
    names(comp_list) <- ct_names
    out_file <- do.call(rbind, comp_list)
    out_file$ct <- gsub('\\..*', '', rownames(out_file))
    write.csv(out_file, paste0(plot_path, "overlap_num_df.csv"))
    for (ct_id in names(comp_list)) {
      ct_df <- comp_list[[ct_id]]
      colnames(ct_df) <- c("comparison", "genes_num")
      ct_df <- separate(ct_df, comparison, into = c("ref_cond", "other_cond"), remove = F, sep = " - ")
      ct_df$ref_cond <- factor(ct_df$ref_cond, groups_ordered[which(groups_ordered %in% unique(ct_df$ref_cond))])
      ct_df$other_cond <- factor(ct_df$other_cond, groups_ordered[which(groups_ordered %in% unique(ct_df$other_cond))])
      cond_palette <- hue_pal()(length(levels(ct_df$ref_cond)))
      names(cond_palette) <- levels(ct_df$ref_cond)
      pdf(paste0(plot_path, ct_id, "_", sex_id, ".pdf"), width = 15)
      print(
        ggplot(ct_df, aes(other_cond, genes_num, fill=other_cond)) +
          geom_bar(stat = "identity", color="black") +
          facet_wrap(~ref_cond, scales = "free") +
          labs(x="", y="Number of shared genes", fill="Datasets", title = paste(ct_id, sex_id, sep = " - ")) +
          fill_palette(cond_palette) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                panel.spacing.x=unit(0, "lines"),
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
                axis.text.x = element_blank(),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    }
  }
}

# 29. Plot the num of shared genes between one condition and all others
# Input: main directory where to save the plots, list of presence dfs, one per each ct, 
    # the order in which to plot the conditions, and the minimum number of conditions ot have each ref gene
# Return: nothing, saves plot and CSV instead

PlotDEGsOverlapHmp <- function(main_dir, ct_df_list, groups_ordered, min_num_conds=2) {
  `%!in%` <- Negate(`%in%`)
  plot_path <- paste0(main_dir, "DEGs_Overlap_across_conditions/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  filt_df <- do.call(rbind, ct_df_list)
  filt_df$presence <- ifelse(filt_df$presence=="yes", 1, 0)
  filt_df$ct <- gsub('\\..*', '', rownames(filt_df))
  for (sex_id in c("F", "M")) {
    print(sex_id)
    sex_df <- subset(filt_df, sex==sex_id & presence==1)
    ct_vec <- vector()
    comp_list <- list()
    ct_names <- vector()
    for (ct in unique(sex_df$ct)) {
      if (length(unique(sex_df[which(sex_df$ct==ct), "condition"]))>2) {
        ct_names <- c(ct_names, ct)
        print(ct)
        ct_cond <- groups_ordered[which(groups_ordered %in% unique(sex_df[which(sex_df$ct==ct), "condition"]))]
        ct_df <- data.frame()
        for (cond in ct_cond) {
          ref_genes <- sex_df[which(sex_df$ct==ct & sex_df$condition==cond), "gene_id"]
          comp_vec <- vector()
          ref_presence <- vector()
          for (other_cond  in ct_cond[!ct_cond == cond]) {
            ref_presence <- c(ref_presence, ifelse(ref_genes %in% sex_df[which(sex_df$ct==ct & sex_df$condition==other_cond), "gene_id"], "Yes", "No"))
            comp_vec <- c(comp_vec, rep(paste(cond, other_cond, sep = " - "), length(ref_genes)))
          }
          ref_genes_ls <- rep(ref_genes, length(ct_cond[!ct_cond == cond]))
          cond_df <- data.frame(comp_vec, ref_genes_ls, ref_presence)
          gene_counts <- as.data.frame(table(cond_df$ref_genes_ls))
          if (length(as.character(gene_counts[which(gene_counts$Freq < min_num_conds),"Var1"])) != 0 & ct != "Dorsal progenitors") {
            genes_to_be_removed <- as.character(gene_counts[which(gene_counts$Freq < min_num_conds),"Var1"])
            cond_df <- subset(cond_df, ref_genes_ls %!in%  genes_to_be_removed)
          }
          ct_df <- rbind(ct_df, cond_df)
        }
        comp_list <- append(comp_list, list(ct_df))
      }
    }
    names(comp_list) <- ct_names
    out_file <- do.call(rbind, comp_list)
    out_file$ct <- gsub('\\..*', '', rownames(out_file))
    write.csv(out_file, paste0(plot_path, "overlap_presence_df.csv"))
    for (ct_id in names(comp_list)) {
      ct_df <- comp_list[[ct_id]]
      colnames(ct_df) <- c("comparison", "gene_id", "presence")
      ct_df <- separate(ct_df, comparison, into = c("ref_cond", "other_cond"), remove = F, sep = " - ")
      ct_df$ref_cond <- factor(ct_df$ref_cond, groups_ordered[which(groups_ordered %in% unique(ct_df$ref_cond))])
      ct_df$other_cond <- factor(ct_df$other_cond, groups_ordered[which(groups_ordered %in% unique(ct_df$other_cond))])
      ct_df$other_cond_presence <- paste(ct_df$other_cond, ct_df$presence, sep=" ")
      ct_df <- ct_df[order(ct_df$other_cond_presence), ]
      pdf(paste0(plot_path, ct_id, "_", sex_id, "_hmp.pdf"), width = 15, height = 15)
      print(
        ggplot(ct_df, aes(other_cond, gene_id, fill=presence)) +
          geom_tile() +
          scale_fill_manual(values = c("Yes"="#F8766D",
                                       "No"="#00BFC4"),
                            guide = guide_legend(reverse = TRUE)) +
          facet_wrap(~ref_cond, scales = "free") +
          labs(x="", y="Genes", fill="Groups", title = paste(ct_id, sex_id, sep = " - ")) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                panel.spacing.x=unit(0, "lines"),
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    }
  }
}


NormDf <- function(ls_norm, common_annot) {
  norm_ls <- list()
  for (sex in c("F", "M")) {
    sex_ls <- list()
    for (id in names(ls_norm)) {
      proj_sex <- do.call(rbind, ls_norm[[id]][[sex]])
      proj_sex <- cbind("ct" = gsub("\\..*", "", rownames(proj_sex)), proj_sex)
      rownames(proj_sex) <- NULL
      proj_sex$common_annot <- rep(NA, nrow(proj_sex))
      for (ct in unique(proj_sex)) {
        proj_sex[which(proj_sex$ct==ct), "common_annot"] <- common_annot[ct]
      }
      sex_ls <- append(sex_ls, list(proj_sex))
    }
    names(sex_ls) <- names(ls_norm)
    sex_ls <- do.call(rbind, sex_ls)
    sex_ls <- cbind("proj_id" = gsub("\\..*", "", rownames(sex_ls)), sex_ls)
    rownames(sex_ls) <- NULL
    norm_ls <- append(norm_ls, list(sex_ls))
  }
  names(norm_ls) <- c("F", "M")
  return(norm_ls)
}


CalcCommonGenes <- function(norm_ls, ct_col) {
  ct_df <- data.frame()
  for (sex in c("F", "M")) {
    for (ct in unique(norm_ls[[sex]][, ct_col])) {
      projs_ct <- unique(norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct), "proj_id"])
      if (length(projs_ct) == 2) {
        ct_df <- rbind(ct_df, c(ct, sex, length(Reduce(intersect, list(norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct & norm_ls[[sex]]$proj_id==projs_ct[1]), "Genes"], 
                                                                         norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct & norm_ls[[sex]]$proj_id==projs_ct[2]), "Genes"])))))
      } else if (length(projs_ct) == 3) {
        ct_df <- rbind(ct_df, c(ct, sex, length(Reduce(intersect, list(norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct & norm_ls[[sex]]$proj_id==projs_ct[1]), "Genes"], 
                                                                         norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct & norm_ls[[sex]]$proj_id==projs_ct[2]), "Genes"],
                                                                         norm_ls[[sex]][which(norm_ls[[sex]][, ct_col]==ct & norm_ls[[sex]]$proj_id==projs_ct[3]), "Genes"])))))
      } else {
        print(paste0(ct, " present only in one project"))
      }
    }
  }
  colnames(ct_df) <- c("ct", "sex", "num_int_genes")
  ct_df$num_int_genes <- as.numeric(ct_df$num_int_genes)
  return(ct_df)
}

PlotCommonGenes <- function(main_dir, ct_df, folder_name, ct_type) {
  plot_path <- paste0(main_dir, folder_name, "_common_genes/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  if (ct_type=="both") {
    pdf(paste0(plot_path, "annot_faceted.pdf"))
    print(
      ggplot(ct_df, aes(ct, num_int_genes, fill=sex)) +
        geom_bar(stat="identity", color="black", position = "dodge") +
        labs(x="Cell types", y="Number of common genes", fill="Sex") +
        scale_y_continuous(breaks=seq(min(ct_df$num_int_genes), max(ct_df$num_int_genes), by=2)) +
        facet_wrap(~annot_type, nrow = 2, scales = "free") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              panel.spacing.x=unit(0, "lines"),
              plot.title = element_text(size=12, face="bold", colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              axis.text.y = element_text(size=8, colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
      
    )
    dev.off()
  } else {
    pdf(paste0(plot_path, ct_type, ".pdf"))
    print(
      ggplot(ct_df, aes(ct, num_int_genes, fill=sex)) +
        geom_bar(stat="identity", color="black", position = "dodge") +
        labs(x=paste0(ct_type, " cell types"), y="Number of common genes", fill="Sex") +
        scale_y_continuous(breaks=seq(min(ct_df$num_int_genes), max(ct_df$num_int_genes), by=2)) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              panel.spacing.x=unit(0, "lines"),
              plot.title = element_text(size=12, face="bold", colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              axis.text.y = element_text(size=8, colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
      
    )
    dev.off()
  }
  
}
