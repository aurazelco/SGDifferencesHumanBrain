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
#install.packages("ggVennDiagram")

library(readxl) # to import xlsx files
library(stringr) # to modify and harmonize names
library(dplyr) # to extract the column from the O'Brien reference
library(RColorBrewer) # to import the Set2 palette
library(VennDiagram) # to plot the Venn Diagram
library(ggplot2) # to plot the percentages of overlap
library(ggVennDiagram) # to plot the second part of Venn diagrams

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

# 2. Function to filter out ns genes and too low FC, and order based on FC 
  # Input: dataframe of DEGs
  # Return: gene list of significant genes as data.frame

Filter_gene <- function( order.gene.df, pval,FC, min_genes) {
  logFC <- log2(FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val_adj"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  #If there are sig genes add index number for each sig gene
  if(nrow(gene.sig) > min_genes) {
    gene.sig$index <- seq.int(nrow(gene.sig))
  } else {
    gene.sig <- data.frame()
  }
  return(gene.sig)
}

# 3. Import All DEGs from F and M for all ct
  # Input: directory where to find ct sub-folders, list of projects ids, the individual project id to look for, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, pval, FC, min_genes=10, ext, row_col) {
  path <- paste0(main_dir, "/01A_DEGs")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    deg <- lapply(deg, function(x) Filter_gene(x, pval, FC, min_genes))
    deg <- deg[lapply(deg,length)>0]
    if (length(deg)==2) {
      for (i in names(deg)) {
        deg_ct <- data.frame("Genes" = rownames(deg[[i]]))
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

ImportDataset <- function(main_dir, folder_list, pval, FC, min_genes=10) {
  ds_list <- list()
  group_names <- vector()
  for (folder in folder_list) {
    ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder, "/outputs"), pval = pval, FC = FC, min_genes=min_genes)))
    group_names <- c(group_names, folder)
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

# 4. Plot the Venn Diagram for each sex
  # Input: main directory, sex, list of genes to compare, the sex, and the palette to be used
  # Return: nothing, plots are saved instead

PlotVennSex <- function(plot_dir, sex_list, sex, venn_col) {
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

# 6. Plots the comparison with a reference df
  # Input: main directory, the reference df, the column in the reference df where the gene names are, 
    # the list of degs imported, the reference name, and parts of string to be removed
  # Return: nothing, plots and CSVs are saved instead

Venn2ndTrimSex <- function(main_dir, ref_df, n_col=2, list_degs, ref_name, to_remove="empty") {
  if (length(list_degs)>4) {
    print("The Venn diagram can only be used with max of 5 groups")
  } else {
    venn_col <- brewer.pal(length(list_degs)+1, "Set2")
    plot_path <- paste0(main_dir, "VennDiagram_2nd_trimester/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    overlap <- list()
    for (sex in c("F", "M")) {
      sex_list <- list(pull(ref_df, n_col))
      for (ds in names(list_degs)) {
        ds_genes <- do.call(rbind,list_degs[[ds]][[sex]])
        #ds_genes$ct <- gsub('(.*).\\w+', '\\1', rownames(ds_genes))
        ds_genes <- ds_genes[!duplicated(ds_genes), ]
        sex_list <- append(sex_list, list(ds_genes))
      }
      names(sex_list) <- c(ref_name, str_remove_all(names(list_degs), to_remove))
      PlotVennSex(plot_path, sex_list, sex, venn_col)
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
    write.csv(data.frame("comparison"=comp, "sex"= sexes, "gene_id"= genes), paste0(plot_path, "sex_intersected_genes.csv"))
  }
}

# 7. Create a list wth all the groups to be compared
  # Input: the imported list of DEGs, the reference df, the column in the ref_df where to find the gene ids, if strings have to be removed from the name
  # Return: list of genes

CreateList <- function(list_degs, ref_df, ref_name, n_col=2, to_remove="empty") {
  groups <- vector()
  ds_sex_genes <- list(pull(ref_df, n_col))
  for (ds in names(list_degs)) {
    for (sex in names(list_degs[[ds]])) {
      ds_sex_genes <- append(ds_sex_genes, list(unique(unlist(list_degs[[ds]][[sex]]))))
      groups <- c(groups, paste(ds, sex, sep = "_"))
    }
  }
  if (to_remove=="empty") {
    names(ds_sex_genes) <- c(ref_name, groups)
  } else {
    names(ds_sex_genes) <- str_remove_all(c(ref_name, groups), to_remove)
  }
  return(ds_sex_genes)
}

# 8. Plots the comparison with a reference df
  # Input: main directory, the reference df, the list of degs imported, the reference name, 
    #  the groups to be compared, the plot tile, the column in the reference df where the gene names are, 
    # and parts of string to be removed
  # Return: nothing, plots and CSVs are saved instead

Venn2ndTrim <- function(main_dir, ref_df, list_degs, ref_name, groups_to_compare, plot_title, n_col=2, to_remove="empty", cat_names) {
  if (length(groups_to_compare)>5) {
    print("The Venn diagram can only be used with max of 5 groups")
  } else {
    plot_path <- paste0(main_dir, "VennDiagram_2nd_trimester/")
    dir.create(plot_path, showWarnings = F, recursive = T)
    ds_sex_genes <- CreateList(list_degs, ref_df, ref_name, n_col, to_remove)
    comp_ls <- ds_sex_genes[groups_to_compare]
    pdf(paste0(plot_path, plot_title, "_Venn.pdf"), width = 10)
    print(
      ggVennDiagram(comp_ls, label = "both", label_alpha = 100, label_percent_digit = 2, category.names = cat_names) + 
        scale_x_continuous(expand = expansion(mult = .2)) +
        theme(legend.position = "bottom")
    )
    dev.off()
    genes <- c(Reduce(intersect, comp_ls))
    comp <- c(rep("all", length(genes)))
    if (length(groups_to_compare)>2) {
      for (i in 1:(length(comp_ls)-1)) {
        comp_int <- comp_ls[(i+1):length(names(comp_ls))]
        for (k in 1:length(comp_int)) {
          genes <- c(genes, intersect(comp_ls[[i]], comp_int[[k]]))
          comp <- c(comp, rep(paste(names(comp_ls)[i], names(comp_int)[k], sep = " - "), length(intersect(comp_ls[[i]], comp_int[[k]]))))
        }
      }
    }
    if (length(groups_to_compare)==5) {
      for (sex in c("F", "M")) {
        sex_comps <- c(ref_name, groups_to_compare[which(grepl(sex, groups_to_compare))])
        sex_ls <- comp_ls[sex_comps]
        if (length(Reduce(intersect, sex_ls))>0) {
          genes <- c(genes, Reduce(intersect, sex_ls))
          comp <- c(comp, rep(paste(plot_title, sex, sep = " - "), length(Reduce(intersect, sex_ls))))
        }
      }
    }
    write.csv(data.frame("comparison"=comp, "gene_id"= genes), paste0(plot_path, plot_title, "_intersected_genes.csv"))
  }
}
