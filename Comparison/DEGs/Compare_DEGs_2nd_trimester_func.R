# Author: Aura Zelco
# Brief description:
  # This script is used to plot the DEGs from UCSC 2nd trimester and O'Brien as Venn diagrams
# Brief procedure:
  # 1. Reads all DEG csv files from all the different datasets
  # 2. Imports the reference dataset - in this case O'Brien et al. 2019
  # 3. Plots the resulting Venn diagrams, one for each sex
# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # df: dataframe
  # ds: dataset

# OBS: this script is sourced in Compare_DEGs_2nd_trimester.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries
#install.packages("readxl")
#install.packages("stringr")
#install.packages("dplyr")
#install.packages("RColorBrewer")
#install.packages("VennDiagram")

library(readxl) # to import xlsx files
library(stringr) # to modify and harmonize names
library(dplyr) # to extract the column from the O'Brien reference
library(RColorBrewer) # to import the Set2 palette
library(VennDiagram) # to plot the Venn Diagram
library(ggplot2) # to plot the percentages of overlap

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

# 2. Import All DEGs from F and M for all ct
  # Input: directory where to find ct sub-folders, list of projects ids, the individual project id to look for, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, individual_projs=vector(), single_proj="", ext, row_col) {
  path <- paste0(main_dir, "/01B_num_DEGs")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    if (length(individual_projs)>0) {
      if (length(grep(single_proj, names(deg))) > 0) {
        deg <- deg[grepl(single_proj, names(deg),fixed = T)]
      } else {
        next
      }
    }
    for (i in names(deg)) {
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
  names(df_F) <- tolower(names_F)
  names(df_M) <- tolower(names_M)
  df_F <- df_F[lapply(df_F,length)>0]
  df_M <- df_M[lapply(df_M,length)>0]
  if (length(df_F) != 0 & length(df_F) != 0) {
    return(list("F"=df_F, "M"=df_M))
  } else {
    return("empty")
  }
}

# 3. Imports UCSC datasets
  # Input: main directory, sub-folders list, list of projects ids
  # Return: list of condition lists, each containing ct lists divided in F and M

ImportDataset <- function(main_dir, folder_list, individual_projs=vector()) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    if (length(individual_projs)==0) {
      ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder, "/outputs"))))
      group_names <- c(group_names, folder)
    } else {
      for (single_proj in individual_projs) {
        single_proj_list <- list(ImportCt(paste0(main_dir, folder), individual_projs = individual_projs, single_proj = single_proj))
        if (single_proj_list!="empty") {
          ds_list <- append(ds_list, single_proj_list)
          group_names <- c(group_names, paste(folder, single_proj, sep = "_"))
        }
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

# 4. Plot the Venn Diagram for each sex
  # Input: main directory, sex, sex
  # Return: nothing, plots are saved instead

PlotVenn <- function(plot_dir, sex_list, sex, venn_col) {
  venn.diagram(
    # General
    x=sex_list, 
    category.names = names(sex_list),
    filename = paste0(plot_dir, sex, "_Venn.png"),
    disable.logging = T, 
    output = T, 
    na="remove",
    # Title
    main = paste0("2nd trimester - ", sex),
    main.fontfamily  = "sans",
    main.fontface =  "bold",
    main.cex = 0.7,
    # Output features
    imagetype="png" ,
    height = 900, 
    width = 900, 
    resolution = 300,
    # Circles
    lwd = 2,
    lty = 'blank',
    col="black",
    fill = venn_col,
    # Numbers
    cex = .6,
    fontfamily = "sans",
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.fontfamily = "sans",
    rotation = 1
  )
}

# 5. Plot the Venn Diagram for each sex
# Input: main directory, sex, sex
# Return: nothing, plots are saved instead

PlotPercentOverlap <- function(plot_dir, overlap) {
  overlapping_genes <- vector()
  for (sex_id in names(overlap)) {
    for (ds_id in names(overlap[[sex_id]])[-1]) {
      overlapping_genes <- c(overlapping_genes, (length(intersect(overlap[[sex_id]][["O'Brien"]], overlap[[sex_id]][[ds_id]]))) * 100 / length(overlap[[sex_id]][[ds_id]]))
    }
  }
  non_overlapping_genes <- 100 - overlapping_genes
  overlap_df <- data.frame("ds"=rep(names(overlap[[sex_id]])[-1], 4),
                           "sex" = rep(c("F", "M"), each=2, times=2),
                           "overlap" = rep(c("Overlapping", "Non overlapping"), each = length(c(overlapping_genes))),
                           "perc" = c(overlapping_genes, non_overlapping_genes))
  pdf(paste0(plot_dir, "overlap_percentage.pdf"))
  print(
    ggplot(overlap_df, aes(ds, perc, fill=sex)) +
      geom_bar(stat="identity", position = "dodge") +
      facet_wrap(~overlap, scales = "free") +
      labs(x="2nd trimester datasets", y="Percentage of DEGs", fill="Sexes") +
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
  )
  dev.off()
}  

# 6. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, the reference df, the list of degs imported
  # Return: nothing, plots are saved instead

Venn2ndTrim <- function(main_dir, ref_df, n_col=2, list_degs, ref_name, to_remove="empty") {
  if (length(list_degs)>4) {
    print("The Venn diagram can only be used with max of 5 groups")
  } else {
    venn_col <- brewer.pal(length(list_degs)+1, "Set2")
    plot_path <- paste0(main_dir, "VennDiagram_2nd_trimester/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    overlap <- list()
    for (sex in c("F", "M")) {
      sex_list <- list(pull(ref_df, 2))
      for (ds in names(list_degs)) {
        ds_genes <- do.call(rbind,list_degs[[ds]][[sex]])
        #ds_genes$ct <- gsub('(.*).\\w+', '\\1', rownames(ds_genes))
        ds_genes <- ds_genes[!duplicated(ds_genes), ]
        sex_list <- append(sex_list, list(ds_genes))
      }
      names(sex_list) <- c(ref_name, str_remove_all(names(list_degs), to_remove))
      #PlotVenn(plot_path, sex_list, sex, venn_col)
      overlap <- append(overlap, list(sex_list))
    }
    names(overlap) <- c("F", "M")
    PlotPercentOverlap(plot_path, overlap)
    genes <- vector()
    comp <- vector()
    sexes <- vector()
    for (sex_id in names(overlap)) {
      genes <- c(genes, Reduce(intersect, overlap[[sex_id]]))
      comp <- c(comp, rep(paste("all", sex_id, sep = "_"), length(Reduce(intersect, overlap[[sex_id]]))))
      sexes <- c(sexes, rep(sex_id, length(Reduce(intersect, overlap[[sex_id]]))))
      for (i in names(overlap[[sex_id]])[-1]) {
        genes <- c(genes, intersect(overlap[[sex_id]][[1]], overlap[[sex_id]][[i]]))
        comp <- c(comp, rep(paste(names(overlap[[sex_id]])[1], i, sep = " - "), 
                            length(intersect(overlap[[sex_id]][[1]], overlap[[sex_id]][[i]]))))
        sexes <- c(sexes, rep(sex_id, 
                            length(intersect(overlap[[sex_id]][[1]], overlap[[sex_id]][[i]]))))
      }
    }
    write.csv(data.frame("comparison"=comp, "sex"= sexes, "gene_id"= genes), paste0(plot_path, "Intersected_genes.csv"))
  }
}

