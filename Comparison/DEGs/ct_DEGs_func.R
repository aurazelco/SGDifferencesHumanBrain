# Author: Aura Zelco
# Brief description:
    # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
    # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
    # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
    # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
    # 4. Plots how many genes are found in all age groups, in all but one, etc
# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # ct: celltype
  # df: dataframe
  # ds: dataset
  # hmps: heatmaps

# OBS: this script is sourced in ct_DEGs.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
library(ggpubr) # to assemble plots together before saving

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
  # Input: directory where to find ct sub-folders, if UCSC or not, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, UCSC_flag="no", ext, row_col) {
  path <- paste0(main_dir, "/01B_num_DEGs")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    for (i in names(deg)) {
      if (UCSC_flag=="no") {
        colnames(deg[[i]]) <- c("Genes")
        rownames(deg[[i]]) <- NULL
        if (grepl("F", i, fixed=TRUE)){
          df_F <- append(df_F, list(deg[[i]]))
          names_F <- c(names_F, sub_ct[ct])
        } else {
          df_M <- append(df_M, list(deg[[i]]))
          names_M <- c(names_M, sub_ct[ct])
        }
      } else {
        deg_ct <- as.data.frame(rownames(deg[[i]]))
        colnames(deg_ct) <- c("Genes")
        rownames(deg_ct) <- NULL
        if (grepl("F", i, fixed=TRUE)){
          df_F <- append(df_F, list(deg_ct))
          names_F <- c(names_F, sub_ct[ct])
        } else {
          df_M <- append(df_M, list(deg_ct))
          names_M <- c(names_M, sub_ct[ct])
        }
      }
    }
  }
  names(df_F) <- tolower(names_F)
  names(df_M) <- tolower(names_M)
  return(list("F"=df_F, "M"=df_M))
}

# 3. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not
  # Return: list of age/condition lists, each containing ct lists divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no") {
  ds_list <- list()
  ct_list <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder))))
      ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
    } else {
      ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder, "/outputs"), UCSC_flag="yes")))
      ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/outputs/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
    }
  }
  names(ds_list) <- folder_list
  return(list("genes"=ds_list, "ct"=unique(ct_list)))
}

# 4. Creates the df for the input ct so that we know if a DEG is found in a certain age or not -> used to generate hmps
  # Input: list of ct dfs, which sex and ct to analyze
  # Return: df with info whether each gene is present in all age/condition groups in which the ct is found

CreatePresenceCtDf <- function(sex_dfs, sex, ct) {
  sub_ct <- sex_dfs[[sex]][which(sex_dfs[[sex]]$common_annot==ct), ]
  ct_sex <-(rep(unique(sub_ct$gene_id), length(unique(sub_ct$age))))
  ct_sex <- cbind(as.data.frame(ct_sex), rep(unique(sub_ct$age), each=length(unique(sub_ct$gene_id))))
  colnames(ct_sex) <- c("gene_id", "age")
  ct_sex$sex <- rep(sex, nrow(ct_sex))
  ct_sex$presence <- rep("no", nrow(ct_sex))
  for (id in unique(ct_sex$gene_id)) {
    for (age_id in unique(ct_sex$age)) {
      if (id %in% sub_ct[which(sub_ct$age==age_id), "gene_id"]) {
        ct_sex[which(ct_sex$age==age_id & ct_sex$gene_id==id), "presence"] <- "yes"
      }
    }
  }
  return(ct_sex)
}

# 5. Creates all PresenceDfs for all cts
  # Input: list of lists, each list corresponding to a specific age/condition-sex-ct combination
  # Return: list of ct dfs, with information on presence of each gene across all ages/conditions

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
        sex_df <- separate(sex_df, groups, into=c("age", "sex", "ct", "gene_num"), sep="\\.")
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

# Plots the heatmap for a certain ct and sex combo
  # Input: presence df, which sex to plot, and the order in which to plot the ages/conditions
  # Return: heatmap plot

PlotDEGsAges <- function(ct_df, sex, age_ordered) {
  ct_df <- ct_df[which(ct_df$sex==sex),]
  ct_df_ordered <- age_ordered[which(age_ordered %in% unique(ct_df$age))]
  ct_df$age <- factor(ct_df$age, ct_df_ordered)
  ct_df <- ct_df[order(ct_df$age), ]
  ct_plot <- ggplot(ct_df, aes(age, gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("yes"="#F8766D",
                                 "no"="#00BFC4"),
                      guide = guide_legend(reverse = TRUE)) +
    labs(x="Developmental Ages", y="Genes", fill="Genes found", title = sex) +
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

# Generates the presence hmps for each ct, putting together F and M from the same ct and saving it as a pdf
  # Input: main directory where to save the plots, the list of presence dfs, and the order in which to plot the ages/conditions
  # Return: nothing, saves the plot instead

PlotCts <- function(main_dir, ct_df_list, age_ordered) {
  plot_path <- paste0(main_dir, "Hmp_DEGs_across_ages/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct in names(ct_df_list)) {
    print(ct)
    f_plot <- PlotDEGsAges(ct_df_list[[ct]], "F", age_ordered)
    m_plot <- PlotDEGsAges(ct_df_list[[ct]], "M", age_ordered)
    ct_plot <- ggarrange(f_plot, m_plot, common.legend = T, legend = "bottom")
    pdf(paste0(plot_path, ct, ".pdf"))
    print(ct_plot)
    dev.off()
  }
}

# Creates dfs which counts in how many ages we find each gene, per sex and ct combo
  # Input: ct df obtained previously, and the sex to analyze
  # Return: df containing for each gene the number of ages which had that gene in their DEGs

GroupsSharingGenes <- function(ct_df, sex_id) {
  ct_df <- ct_df[which(ct_df$sex==sex_id),]
  gene_id <- vector()
  age_count <- vector()
  for (id_gene in unique(ct_df$gene_id)) {
    gene_id <- c(gene_id, id_gene)
    age_count <- c(age_count, length(unique(ct_df[which(ct_df$gene_id==id_gene & ct_df$presence=="yes"), "age"])))
  }
  sex <- rep(sex_id, length(gene_id))
  return(data.frame(gene_id, sex, age_count))
}

# Create Count Dfs for all cts
  # Input: presence df list
  # Return: list of df containing the number of ages for each gene, for each ct

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

# Plot count dfs for each ct
  # Input: the count df of one ct
  # Return: bar plot of how many genes are shared among how many groups

PlotCountDfCt <- function(ct_df) {
  ct_plot <- ggplot(ct_df, aes(age_count, fill=as.factor(age_count))) +
    geom_bar() +
    facet_wrap(~sex, scales = "free") +
    scale_x_continuous(breaks=seq(min(ct_df$age_count), max(ct_df$age_count),by=1)) +
    labs(y="Gene absolute count", fill="Number of ages sharing genes") +
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

# Plot all count dfs for all cts
  # Input: main directory where to save the plots, the list of count dfs
  # Return: nothing, saves the plot instead

PlotCountCt <- function(main_dir, gene_count_dfs) {
  plot_path <- paste0(main_dir, "Num_DEGs_across_ages/")
  dir.create(plot_path, showWarnings = F, recursive = T)
  for (ct in names(gene_count_dfs)) {
    print(ct)
    ct_plot <- PlotCountDfCt(gene_count_dfs[[ct]])
    pdf(paste0(plot_path, ct, ".pdf"))
    print(ct_plot)
    dev.off()
  }
}