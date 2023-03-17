# Author: Aura Zelco
# Brief description:
  # This script is used to run the enrichment analysis on multiple databases/terms, all using the DEGs from the projects
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all groups
  # 3. For each enrichments, compares F v M in each ct-groups combo and in each sex-ct combo across groups
    # clusterProfile:
      # a. GO - Gene Ontology 
      # b. KEGG - Kyoto Encyclopedia of Genes and Genomes
      # c. DO - Disease Ontology
      # d. DGN - DisGeNET
    # enrichR:
      # a. DSigDB
      # b. GWAS_Catalog_2019
    # disgenet2r:
      # a. Curated DisGeNET
  # 4. Saves the plots and CSV results

# Documentation abbreviations:
  # deg: differentially expressed genes
  # F and M: females and males
  # ct: celltype
  # df: dataframe


# OBS: this script is sourced in Compare_Enrichment.R

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
# BiocManager::install("DOSE")
# install.packages("tidyr")
# install.packages("ggnewscale")
# install.packages("ggplot2")
# install.packages("stringr")
# install.packages("enrichR")
# install.packages("devtools")
# library(devtools)
# How to install SPARQL from source since it has been archived from CRAN
# url <- "https://cran.r-project.org/src/contrib/Archive/SPARQL/SPARQL_1.16.tar.gz"
# pkgFile <- "SPARQL_1.16.tar.gz"
# R_libraries <- "/Library/Frameworks/R.framework/Versions/4.1-arm64/Resources/library/"
# download.file(url = url, destfile = paste0(R_libraries, pkgFile))
# Expand the zip file using whatever system functions are preferred
# $ tar -xf SPARQL_1.16.tar.g
# look at the DESCRIPTION file in the expanded package directory
# Install dependencies list in the DESCRIPTION file
# install.packages(c("XML", "RCurl"))
# Install package
# install.packages(pkgs=paste0(R_libraries, pkgFile), type="source", repos=NULL)
# Delete package tarball
# unlink(pkgFile)
# install_bitbucket("ibi_group/disgenet2r")
# install.packages("readxl")

library(clusterProfiler) # to run the enrichment analysis (GO, KEGG)
library(DOSE)  # to run the enrichment analysis (DO)
library(org.Hs.eg.db) # the organism database to use
library(tidyr) # to re-arrange dfs
library(ggnewscale) # to plot the cnet plots
library(tidyr) # to re-arrange dfs
library(ggplot2) # to plot
library(stringr) # to format strings
library(enrichR) # database
library(disgenet2r) # database
library(readxl) # to read excel files
library(dplyr) # to re-arrange dfs
library(biomaRt) # to query to which chromosome the shared terms belong to


# Sets the EnrichR database for Human genes
setEnrichrSite("Enrichr") 

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

# 5. Creates a list, each element the gene list from each groups present for a specific ct and sex
  # Input: the dataframe, the ct and sex to be analyzed, the order in which the groups should be plotted
  # Return: the list of gene lists, from a specific ct and sex, for each groups

ExtractSexCt <- function(sex_df, ct, sex, groups_ordered) {
  sex_ct <- sex_df[which(sex_df$common_annot==ct & sex_df$sex==sex), ]
  sex_ct <- split(sex_ct$gene_id, sex_ct$groups)
  sex_order <- groups_ordered[which(groups_ordered %in% names(sex_ct))]
  sex_ct <- sex_ct[sex_order]
  return(sex_ct)
}

# 6. Creates a list, each element the gene list from each cts present for a specific group and sex
  # Input: the dataframe, the group and sex to be analyzed, the order in which the cts should be plotted
  # Return: the list of gene lists, from a specific group and sex, for each ct

ExtractSexGroup <- function(sex_df, group_id, sex, cts_ordered) {
  sex_gr <- sex_df[which(sex_df$groups==group_id & sex_df$sex==sex), ]
  sex_gr <- split(sex_gr$gene_id, sex_gr$common_annot)
  sex_order <- cts_ordered[which(cts_ordered %in% names(sex_gr))]
  sex_gr <- sex_gr[sex_order]
  return(sex_gr)
}

# 7. Compares the GOs of a list of genes
  # Input: list of genes to be compared, which GO (BP, MO or CC),  
    # if a minimum threshold for how many genes in each module should be used,
  # Return: enriched GO (formal class compareClusterResult)

compareGO <- function(sex_list, GO_ont, gene_thresh="no"){
  enrich <- compareCluster(geneCluster  = sex_list, 
                           fun          = "enrichGO",
                           keyType      = "SYMBOL",
                           ont          = GO_ont,
                           OrgDb        = "org.Hs.eg.db",
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  if (is.numeric(gene_thresh)) {
    enrich <- gsfilter(enrich, min=gene_thresh)
  }
  return(enrich)
}

# 8. finds the EntrezID for a gene
  # Input: gene symbol
  # Return: EntrezID(s) for the gene

GenetoENTREZ <- function(symbol){
  library(org.Hs.eg.db)
  hs <- org.Hs.eg.db
  entrez.df <- AnnotationDbi::select(hs, 
                                     keys = symbol,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")
  return(entrez.df)
}

# 9. Compares the KEGG pathways of a list of genes
  # Input: list of genes to be compared, if a minimum threshold for how many genes in each module should be used
  # Return: enriched KEGG (formal class compareClusterResult)

compareKEGG <- function(sex_list, gene_thresh="no"){
  sex_list_kegg <- lapply(sex_list, function(x) 
                          {gene.df <- GenetoENTREZ(x)
                          gene.df <- gene.df[which(!is.na(gene.df$ENTREZID)), ]
                          return(gene.df$ENTREZID)})
  names(sex_list_kegg) <- names(sex_list)
  enrich <- compareCluster(geneCluster  = sex_list_kegg, 
                           fun          = "enrichKEGG",
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05
                           )
  if (is.numeric(gene_thresh)) {
    enrich <- gsfilter(enrich, min=gene_thresh)
  }
  return(enrich)
}

# 10. Compares the DO of a list of genes
  # Input: list of genes to be compared, if a minimum threshold for how many genes in each module should be used
  # Return: enriched KEGG (formal class compareClusterResult)

compareDO <- function(sex_list, gene_thresh="no"){
  sex_list_DO <- lapply(sex_list, function(x) 
                        {gene.df <- GenetoENTREZ(x)
                        gene.df <- gene.df[which(!is.na(gene.df$ENTREZID)), ]
                        return(gene.df$ENTREZID)})
  names(sex_list_DO) <- names(sex_list)
  enrich <- compareCluster(geneCluster  = sex_list_DO, 
                           fun          = "enrichDO",
                           ont          = "DO",
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
  if (is.numeric(gene_thresh)) {
    enrich <- gsfilter(enrich, min=gene_thresh)
  }
  return(enrich)
}

# 11. Compares the DGN of a list of genes
  # Input: list of genes to be compared, if a minimum threshold for how many genes in each module should be used
  # Return: enriched KEGG (formal class compareClusterResult)

compareDGN <- function(sex_list, gene_thresh="no"){
  sex_list_DGN <- lapply(sex_list, function(x) 
                        {gene.df <- GenetoENTREZ(x)
                        gene.df <- gene.df[which(!is.na(gene.df$ENTREZID)), ]
                        return(gene.df$ENTREZID)})
  names(sex_list_DGN) <- names(sex_list)
  enrich <- compareCluster(geneCluster  = sex_list_DGN, 
                            fun          = "enrichDGN",
                            pAdjustMethod= "BH",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.05)
  if (is.numeric(gene_thresh)) {
    enrich <- gsfilter(enrich, min=gene_thresh)
  }
  return(enrich)
}

# 12. Compares the selected enrichment module of a list of genes
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the enrichment to be analyzed, 
    # the GO to analyze (if GO is the module), and  if a minimum threshold for how many genes in each module should be used
  # Return: nothing, saves plots and CSVs instead

EnrichFvM <- function(main_dir, sex_df, enrich_module, GO_ont="BP", gene_thresh="no") {
  out_path <- paste0(main_comparison, enrich_module, "_comparison_ct_sex/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (cond in unique(sex_df$groups)) {
    print(cond)
    cond_path <- paste0(out_path, cond, "/")
    dir.create(cond_path, recursive = T, showWarnings = F)
    cond_df <- sex_df[which(sex_df$groups==cond), ]
    for (ct in unique(cond_df$common_annot)) {
      print(ct)
      f_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="F"), "gene_id"]
      m_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="M"), "gene_id"]
      sex_comp <- list("F"=f_genes, "M"=m_genes)
      if (enrich_module=="GO") {
        try({
          print(paste0("Calculating the ", GO_ont, " results for ", ct))
          sex_ct_enriched <-  compareGO(sex_comp, GO_ont, gene_thresh)
          png(paste0(cond_path, "GO_", GO_ont, "_", ct, ".png"))
          print(dotplot(sex_ct_enriched, by = 'count', title=ct))
          dev.off()
          csv_ct_enriched <- as.data.frame(sex_ct_enriched)
          print(paste0("Saving the CSV GO results for ", ct))
          write.csv(csv_ct_enriched, paste0(cond_path, "GO_", GO_ont, "_", ct, ".csv"))
        })
      } else if (enrich_module=="KEGG") {
        try({
          print(paste0("Calculating the KEGG results for ", ct))
          sex_ct_enriched <-  compareKEGG(sex_comp, gene_thresh)
          png(paste0(cond_path, "KEGG_", ct, ".png"))
          print(dotplot(sex_ct_enriched, by = 'count', title=ct))
          dev.off()
          csv_ct_enriched <- as.data.frame(sex_ct_enriched)
          print(paste0("Saving the CSV KEGG results for ", ct))
          write.csv(csv_ct_enriched, paste0(cond_path, "KEGG_", ct, ".csv"))
        })
      } else if (enrich_module=="DO") {
        try({
          print(paste0("Calculating the DO results for ", ct))
          sex_ct_enriched <-  compareDO(sex_comp, gene_thresh)
          png(paste0(cond_path, "DO_", ct, ".png"))
          print(dotplot(sex_ct_enriched, by = 'count', title=ct))
          dev.off()
          csv_ct_enriched <- as.data.frame(sex_ct_enriched)
          print(paste0("Saving the CSV DO results for ", ct))
          write.csv(csv_ct_enriched, paste0(cond_path, "DO_", ct, ".csv"))
        })
      } else if (enrich_module=="DGN") {
        try({
          print(paste0("Calculating the DGN results for ", ct))
          sex_ct_enriched <-  compareDGN(sex_comp, gene_thresh)
          png(paste0(cond_path, "DGN_", ct, ".png"))
          print(dotplot(sex_ct_enriched, by = 'count', title=ct))
          dev.off()
          csv_ct_enriched <- as.data.frame(sex_ct_enriched)
          print(paste0("Saving the CSV DGN results for ", ct))
          write.csv(csv_ct_enriched, paste0(cond_path, "DGN_", ct, ".csv"))
        })
      }
    }
  }
}


# 13. Compares the DEGs from all groups in a specific ct-sex group
  # Input: main directory where to save the plots, the dataframe containing all DEGs, , the enrichment to be analyzed, 
    # the GO to analyze (if GO is the module), if a minimum threshold for how many genes in each module should be used, 
    # the order in which plot the groups, and if the x-axis labels should be rotated by 90 degrees, the threshold for the adjusted p-value to use
  # Return: nothing, saves plots and CSVs instead

EnrichCts <- function(main_dir, sex_df, enrich_module, GO_ont="BP", gene_thresh="no", groups_ordered, rotate_x_axis=F, adj_pval_thresh=0.05) {
  out_path <- paste0(main_comparison, enrich_module, "_comparison_cts/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (ct in unique(sex_df$common_annot)) {
    print(ct)
    ct_path <- paste0(out_path, ct, "/")
    dir.create(ct_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      sex_ct <- ExtractSexCt(sex_df, ct, sex, groups_ordered)
      if (enrich_module=="GO") {
        try({
          print(paste0("Calculating the ", GO_ont, " results for ", sex))
          sex_ct_GO <-  compareGO(sex_ct, GO_ont, gene_thresh)
          csv_ct_GO <- as.data.frame(sex_ct_GO)
          print(paste0("Saving the CSV GO results for ", sex))
          write.csv(csv_ct_GO, paste0(ct_path, "GO_", GO_ont, "_", sex, ".csv"))
          sex_ct_GO@compareClusterResult <- sex_ct_GO@compareClusterResult[which(sex_ct_GO@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_ct_plot <- dotplot(sex_ct_GO, by = 'count', title=paste(ct, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_ct_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "GO_", GO_ont, "_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_ct_plot)
          dev.off()
        })
      } else if (enrich_module=="KEGG") {
        try({
          print(paste0("Calculating the KEGG results for ", sex))
          sex_ct_KEGG <- compareKEGG(sex_ct, gene_thresh)
          return(sex_ct_KEGG)
          csv_ct_KEGG <- as.data.frame(sex_ct_KEGG)
          print(paste0("Saving the CSV KEGG results for ", sex))
          write.csv(csv_ct_KEGG, paste0(ct_path, "KEGG_", sex, ".csv"))
          sex_ct_KEGG@compareClusterResult <- sex_ct_KEGG@compareClusterResult[which(sex_ct_KEGG@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_ct_plot <- dotplot(sex_ct_KEGG, by = 'count', title=paste(ct, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_ct_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "KEGG_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_ct_plot)
          dev.off()
          png(paste0(ct_path, "KEGG_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(setReadable(sex_ct_KEGG, 'org.Hs.eg.db', 'ENTREZID'),cex_label_category = 1.2))
          dev.off()
        })
      } else if (enrich_module=="DO") {
        try({
          print(paste0("Calculating the DO results for ", sex))
          sex_ct_DO <-  compareDO(sex_ct, gene_thresh)
          csv_ct_DO <- as.data.frame(sex_ct_DO)
          print(paste0("Saving the CSV DO results for ", sex))
          write.csv(csv_ct_DO, paste0(ct_path, "DO_", sex, ".csv"))
          sex_ct_DO@compareClusterResult <- sex_ct_DO@compareClusterResult[which(sex_ct_DO@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_ct_plot <- dotplot(sex_ct_DO, by = 'count', title=paste(ct, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_ct_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "DO_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_ct_plot)
          dev.off()
          png(paste0(ct_path, "DO_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(sex_ct_DO, node_label="category",cex_label_category = 1.2))
          dev.off()
        })
      } else if (enrich_module=="DGN") {
        try({
        print(paste0("Calculating the DGN results for ", sex))
        sex_ct_DGN <-  compareDGN(sex_ct, gene_thresh)
        csv_ct_DGN <- as.data.frame(sex_ct_DGN)
        print(paste0("Saving the CSV DGN results for ", sex))
        write.csv(csv_ct_DGN, paste0(ct_path, "DGN_", sex, ".csv"))
        sex_ct_DGN@compareClusterResult <- sex_ct_DGN@compareClusterResult[which(sex_ct_DGN@compareClusterResult$p.adjust<=adj_pval_thresh),]
        sex_ct_plot <- dotplot(sex_ct_DGN, by = 'count', title=paste(ct, sex, sep=" - "), showCategory=2)
        if (rotate_x_axis) {sex_ct_plot$theme$axis.text.x$angle <- 90}
        png(paste0(ct_path, "DGN_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
        print(sex_ct_plot)
        dev.off()
        })
      }
    }
  }
}

# 14. Compares the DEGs from all groups in a specific ct-sex group
  # Input: main directory where to save the plots, the dataframe containing all DEGs, , the enrichment to be analyzed, 
    # the GO to analyze (if GO is the module), if a minimum threshold for how many genes in each module should be used, 
    # the order in which plot the cts, and if the x-axis labels should be rotated by 90 degrees, the threshold for the adjusted p-value to use
  # Return: nothing, saves plots and CSVs instead

EnrichGroups <- function(main_dir, sex_df, enrich_module, GO_ont="BP", gene_thresh="no", cts_ordered, rotate_x_axis=F, adj_pval_thresh=0.05) {
  out_path <- paste0(main_comparison, enrich_module, "_comparison_groups/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (group in unique(sex_df$groups)) {
    print(group)
    group_path <- paste0(out_path, group, "/")
    dir.create(group_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      sex_group <- ExtractSexGroup(sex_df, group, sex, cts_ordered)
      if (enrich_module=="GO") {
        try({
          print(paste0("Calculating the ", GO_ont, " results for ", sex))
          sex_group_GO <-  compareGO(sex_group, GO_ont, gene_thresh)
          csv_group_GO <- as.data.frame(sex_group_GO)
          print(paste0("Saving the CSV GO results for ", sex))
          write.csv(csv_group_GO, paste0(group_path, "GO_", GO_ont, "_", sex, ".csv"))
          sex_group_GO@compareClusterResult <- sex_group_GO@compareClusterResult[which(sex_group_GO@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_group_plot <- dotplot(sex_group_GO, by = 'count', title=paste(group, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_group_plot$theme$axis.text.x$angle <- 90}
          png(paste0(group_path, "GO_", GO_ont, "_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_group_plot)
          dev.off()
        })
      } else if (enrich_module=="KEGG") {
        try({
          print(paste0("Calculating the KEGG results for ", sex))
          sex_group_KEGG <- compareKEGG(sex_group, gene_thresh)
          csv_group_KEGG <- as.data.frame(sex_group_KEGG)
          print(paste0("Saving the CSV KEGG results for ", sex))
          write.csv(csv_group_KEGG, paste0(group_path, "KEGG_", sex, ".csv"))
          sex_group_KEGG@compareClusterResult <- sex_group_KEGG@compareClusterResult[which(sex_group_KEGG@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_group_plot <- dotplot(sex_group_KEGG, by = 'count', title=paste(group, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_group_plot$theme$axis.text.x$angle <- 90}
          png(paste0(group_path, "KEGG_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_group_plot)
          dev.off()
          png(paste0(group_path, "KEGG_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(setReadable(sex_group_KEGG, 'org.Hs.eg.db', 'ENTREZID'),cex_label_category = 1.2))
          dev.off()
        })
      } else if (enrich_module=="DO") {
        try({
          print(paste0("Calculating the DO results for ", sex))
          sex_group_DO <-  compareDO(sex_group, gene_thresh)
          csv_group_DO <- as.data.frame(sex_group_DO)
          print(paste0("Saving the CSV DO results for ", sex))
          write.csv(csv_group_DO, paste0(group_path, "DO_", sex, ".csv"))
          sex_group_DO@compareClusterResult <- sex_group_DO@compareClusterResult[which(sex_group_DO@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_group_plot <- dotplot(sex_group_DO, by = 'count', title=paste(group, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_group_plot$theme$axis.text.x$angle <- 90}
          png(paste0(group_path, "DO_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_group_plot)
          dev.off()
          png(paste0(group_path, "DO_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(sex_group_DO, node_label="category",cex_label_category = 1.2))
          dev.off()
        })
      } else if (enrich_module=="DGN") {
        try({
          print(paste0("Calculating the DGN results for ", sex))
          sex_group_DGN <-  compareDGN(sex_group, gene_thresh)
          csv_group_DGN <- as.data.frame(sex_group_DGN)
          print(paste0("Saving the CSV DGN results for ", sex))
          write.csv(csv_group_DGN, paste0(group_path, "DGN_", sex, ".csv"))
          sex_group_DGN@compareClusterResult <- sex_group_DGN@compareClusterResult[which(sex_group_DGN@compareClusterResult$p.adjust<=adj_pval_thresh),]
          sex_group_plot <- dotplot(sex_group_DGN, by = 'count', title=paste(group, sex, sep=" - "), showCategory=2)
          if (rotate_x_axis) {sex_group_plot$theme$axis.text.x$angle <- 90}
          png(paste0(group_path, "DGN_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_group_plot)
          dev.off()
        })
      }
    }
  }
}


# 15. Searches the selected database in enrichR
  # Input: the gene list, the database to look into
  # Return: dataframe with the retrieved information

EnrichR_fun <- function(gene_ls, dbsx){
  enrichR_obj <- enrichr(gene_ls, dbsx)
  enrichR_df <- enrichR_obj[[1]]
  gene_count <-  unlist(lapply(enrichR_df$Overlap, function(x) { unlist(strsplit(x, "/"))[1] } ))
  enrichR_df$gene_count <- gene_count
  if (dbsx=="TRANSFAC_and_JASPAR_PWMs") {
    enrichR_df <- enrichR_df[which(grepl("human", enrichR_df$Term)), ]
  }
  return(enrichR_df)
}

# 16. Searches in the DisGeNET database using disgenet2r
  # Input: the gene list
  # Return: dataframe with the retrieved information

EnrichDisgenet2r_fun <- function(gene_ls){
  dgn_res <-disease_enrichment( entities = gene_ls, 
                                vocabulary = "HGNC", 
                                database = "CURATED")
  dgn_res <- dgn_res@qresult
  return(dgn_res)
}

# 17. Searches in the selected packages and databases for disease-enrichment
  # Input:  list of genes to be compared, the database to look into, the package to use
  # Return: list of enriched terms for each id (group or ct)

EnrichId <- function(sex_id, package, dbsx, name_col){
  if(package == 'EnrichR'){
    enrich_ls <- lapply(sex_id, function(x) EnrichR_fun(x, dbsx))
  }
  else if(package == 'DisGeNET2r'){
    enrich_ls <- lapply(sex_id, function(x) EnrichDisgenet2r_fun(x))
  }
  for (id in names(enrich_ls)) {
    id_df <- enrich_ls[[id]]
    if  (!is.null(dim(id_df)) ){
      id_df[, name_col] <- rep(id, nrow(id_df))
      enrich_ls[[id]] <- id_df
    }
  }
  return(enrich_ls)
}



# 18. Selects the top 5 enriched terms to combine in one df for plotting
  # Input: the enriched list
  # Return: the filtered list of terms

SelectTop5 <- function(enrich_df){
  head_df <- lapply(enrich_df, function(x) head(x, 5))
  combind_df <- Reduce(rbind, head_df)
  # Add gene ratio column
  #combind_df$gene_ratio <- sapply(combind_df$Overlap, function(x) eval(parse(text = x)) )
  return(combind_df)
}

# 19. Calculates the enrichment in the package and database specified for each ct-sex combo across all cts
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
    # the database, the order in which plot the groups
  # Return: nothing, saves plots and CSVs instead

EnrichOtherDB <- function(main_dir, sex_df, package, dbsx, groups_ordered){
  dbsx_path <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  out_path <- paste0(main_dir, package, "_", dbsx_path, "/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (ct in unique(sex_df$common_annot)) {
    print(ct)
    ct_path <- paste0(out_path, ct, "/")
    dir.create(ct_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      filt_flag <- T
      sex_ct <- ExtractSexCt(sex_df, ct, sex, groups_ordered)
      enrich_df <- EnrichId(sex_ct, package, dbsx, 'groups')
      top5 <-SelectTop5(enrich_df)
      write.csv(top5, paste0(ct_path, "top5_", sex,".csv"))
      if (package == 'EnrichR') {
        if (nrow(top5[top5['Adjusted.P.value']< 0.05,] > 1)) {
          filtered_top <- top5[top5['Adjusted.P.value']< 0.05,]
          filtered_top$gene_count <- as.numeric(filtered_top$gene_count)
          if (nrow(filtered_top[filtered_top['gene_count'] > 1,] > 1)) {
            filtered_top <- filtered_top[filtered_top['gene_count'] > 1,]
            y_var = 'Term'
            size_var = 'gene_count'
            color_var = 'Adjusted.P.value'
          } else {
            print("No significant terms with more than 1 gene")
            filt_flag <- F
          }
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      } else if (package == 'DisGeNET2r') {
        if (nrow(top5[top5['FDR'] < 0.05,])) {
          filtered_top <- top5[top5['FDR'] < 0.05,]
          y_var = 'Description'
          size_var = 'Count'
          color_var = 'FDR'
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      }
      if (filt_flag) {
        filtered_top <- filtered_top[, c("groups", y_var, size_var, color_var)]
        colnames(filtered_top) <- c("groups", "term", "gene_count", "adj_pval")
        pdf(paste0(ct_path, sex, ".pdf"))
        print(
          ggplot(filtered_top, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(title = paste(ct, sex,sep = " - "), y = "", x = "Groups", size = "Gene count", color = "Adjusted p-value") +
            scale_size_continuous(range=c(3, 8)) + 
            scale_y_discrete(labels=function(x) str_wrap(x,width=40)) +
            theme(
              plot.title = element_text(size=14, face="bold", colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
              legend.position = "right", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
        )
        dev.off()
      }
    }
  }
}

# 20. Calculates the enrichment in the package and database specified comparing for each ct-groups combo the sexes
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
    # the database, the order in which plot the groups
  # Return: nothing, saves plots and CSVs instead

EnrichOtherDBFvM <- function(main_dir, sex_df, package, dbsx, groups_ordered){
  dbsx_path <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  out_path <- paste0(main_dir, package, "_", dbsx_path, "_ct_sex/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (cond in unique(sex_df$groups)) {
    print(cond)
    cond_df <- sex_df[which(sex_df$groups==cond), ]
    for (ct in unique(cond_df$common_annot)) {
      cond_path <- paste0(out_path, cond, "/", ct, "/")
      dir.create(cond_path, recursive = T, showWarnings = F)
      filt_flag <- T
      print(ct)
      f_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="F"), "gene_id"]
      m_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="M"), "gene_id"]
      sex_comp <- list("F"=f_genes, "M"=m_genes)
      enrich_df <- EnrichId(sex_ct, package, dbsx, 'groups')
      top5 <-SelectTop5(enrich_df)
      write.csv(top5, paste0(cond_path, "top5.csv"))
      if (package == 'EnrichR') {
        if (nrow(top5[top5['Adjusted.P.value']< 0.05,] > 1)) {
          filtered_top <- top5[top5['Adjusted.P.value']< 0.05,]
          filtered_top$gene_count <- as.numeric(filtered_top$gene_count)
          if (nrow(filtered_top[filtered_top['gene_count'] > 1,] > 1)) {
            filtered_top <- filtered_top[filtered_top['gene_count'] > 1,]
            y_var = 'Term'
            size_var = 'gene_count'
            color_var = 'Adjusted.P.value'
          } else {
            print("No significant terms with more than 1 gene")
            filt_flag <- F
          }
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      } else if (package == 'DisGeNET2r') {
        if (nrow(top5[top5['FDR'] < 0.05,])) {
          filtered_top <- top5[top5['FDR'] < 0.05,]
          y_var = 'Description'
          size_var = 'Count'
          color_var = 'FDR'
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      }
      if (filt_flag) {
        filtered_top <- filtered_top[, c("groups", y_var, size_var, color_var)]
        colnames(filtered_top) <- c("groups", "term", "gene_count", "adj_pval")
        pdf(paste0(cond_path, "sexes.pdf"))
        print(
          ggplot(filtered_top, aes(groups, term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(y = "", x = "Sexes", size = "Gene count", color = "Adjusted p-value") +
            scale_size_continuous(range=c(3, 8)) + 
            scale_y_discrete(labels=function(x) str_wrap(x,width=40)) +
            theme(
              plot.title = element_text(size=14, face="bold", colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
              legend.position = "right", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
        )
        dev.off()
      }
    }
  }
}


# 21. Calculates the enrichment in the package and database specified for each ct-sex combo across all cts
# Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
# the database, the order in which plot the groups
# Return: nothing, saves plots and CSVs instead

EnrichOtherDBGroup <- function(main_dir, sex_df, package, dbsx, cts_ordered){
  dbsx_path <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  out_path <- paste0(main_dir, package, "_", dbsx_path, "_groups/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (group in unique(sex_df$groups)) {
    print(group)
    group_path <- paste0(out_path, group, "/")
    dir.create(group_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      filt_flag <- T
      sex_group <- ExtractSexGroup(sex_df, group, sex, cts_ordered)
      enrich_df <- EnrichId(sex_group, package, dbsx, "cts")
      top5 <-SelectTop5(enrich_df)
      write.csv(top5, paste0(group_path, "top5_", sex,".csv"))
      if (package == 'EnrichR') {
        if (nrow(top5[top5['Adjusted.P.value']< 0.05,] > 1)) {
          filtered_top <- top5[top5['Adjusted.P.value']< 0.05,]
          filtered_top$gene_count <- as.numeric(filtered_top$gene_count)
          if (nrow(filtered_top[filtered_top['gene_count'] > 1,] > 1)) {
            filtered_top <- filtered_top[filtered_top['gene_count'] > 1,]
            y_var = 'Term'
            size_var = 'gene_count'
            color_var = 'Adjusted.P.value'
          } else {
            print("No significant terms with more than 1 gene")
            filt_flag <- F
          }
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      } else if (package == 'DisGeNET2r') {
        if (nrow(top5[top5['FDR'] < 0.05,])) {
          filtered_top <- top5[top5['FDR'] < 0.05,]
          y_var = 'Description'
          size_var = 'Count'
          color_var = 'FDR'
        } else {
          print("No significant terms with adjusted p-value < 0.05")
          filt_flag <- F
        }
      }
      if (filt_flag) {
        filtered_top <- filtered_top[, c("cts", y_var, size_var, color_var)]
        colnames(filtered_top) <- c("cts", "term", "gene_count", "adj_pval")
        pdf(paste0(group_path, sex, ".pdf"))
        print(
          ggplot(filtered_top, aes(factor(cts, cts_ordered[which(cts_ordered %in% cts)]), term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(title = paste(group, sex,sep = " - "), y = "", x = "Cell types", size = "Gene count", color = "Adjusted p-value") +
            scale_size_continuous(range=c(3, 8)) + 
            scale_y_discrete(labels=function(x) str_wrap(x,width=40)) +
            theme(
              plot.title = element_text(size=14, face="bold", colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
              legend.position = "right", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
        )
        dev.off()
      }
    }
  }
}


# 22. Plot the presence heatmap
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

# 23. Plot the presence number of genes
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

# 24. Calculates the % of known markers in the DEGs
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

# 25. Plot the presence number of genes as percentage
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

# 26. Plot the presence number of genes as percentage
  # Input: the presence df, the groups order
  # Return: the plot

PlotBarPlotRefPercFaceted <- function(ref_perc, groups_ordered) {
  ref_plot <- ggplot(ref_perc, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), perc, fill = ref_ct)) +
    geom_bar(stat="identity", color="black", position = "dodge") +
    facet_grid(sex ~ ct, scales = "free") +
    labs(x="Groups", y="Markers %", fill="Reference cell types") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black"),
          axis.ticks.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ref_plot)
}


# 27. Plots if gebes from a reference df are found or not in the DEGs
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the reference df, 
    # the order in which plot the groups (and which groups to plot), the vector to use for plot titles
  # Return: nothing, saves plot instead

PlotRefCt <- function(main_dir, sex_df, ref_df, groups_ordered, ref_df_name="ref", plot_titles){
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
  for (ref_ct_id in unique(ref_presence_df$ref_ct)) {
    print(plot_titles[[ref_ct_id]])
    pdf(paste0(out_path, plot_titles[ref_ct_id], "_hmp.pdf"), height = 15, width = 10)
    print(PlotHmpRef(ref_presence_df, ref_ct_id, plot_titles))
    dev.off()
    pdf(paste0(out_path, plot_titles[ref_ct_id], "_barplot.pdf"), height = 4, width = 16)
    print(PlotBarPlotRef(ref_presence_df, ref_ct_id, plot_titles))
    dev.off()
    ref_perc <- RefPerc(ref_presence_df, ref_ct_id, sex_df)
    pdf(paste0(out_path, plot_titles[ref_ct_id], "_barplot_perc.pdf"), height = 4, width = 16)
    print(PlotBarPlotRefPerc(ref_perc, ref_ct_id, plot_titles, groups_ordered))
    dev.off()
  }
  ref_perc <- list()
  ref_perc <- lapply(1:length(unique(ref_presence_df$ref_ct)), function(x) RefPerc(ref_presence_df, unique(ref_presence_df$ref_ct)[x], sex_df, plot_titles))
  ref_perc <- do.call(rbind, ref_perc)
  ref_plot <- PlotBarPlotRefPercFaceted(ref_perc, groups_ordered)
  pdf(paste0(out_path, "Faceted_barplot_perc.pdf"), height = 8, width = 15)
  print(ref_plot)
  dev.off()
}


# 28. Import results from disease-related enrichments
  # Input: main directory where to find the files, which db to use, which comparison to import
  # Return: list of files

ImportDBresults <- function(main_dir, dbsx, which_dbsx) {
  dbsx_path <- paste0(main_dir, dbsx, which_dbsx)
  cts <- list.dirs(dbsx_path, recursive = F, full.names = F)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in cts) {
    dbsx_files <- ImportFiles(paste(dbsx_path, ct, sep = "/"))
    for (i in names(dbsx_files)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(dbsx_files[[i]]))
        names_F <- c(names_F, ct)
      } else {
        df_M <- append(df_M, list(dbsx_files[[i]]))
        names_M <- c(names_M, ct)
      }
    }
  }
  names(df_F) <- names_F
  names(df_M) <- names_M
  dbsx_ls <- list("F"=df_F, "M"=df_M)
  for (sex in names(dbsx_ls)) {
    dbsx_ls[[sex]] <- do.call(rbind, dbsx_ls[[sex]])
    dbsx_ls[[sex]] <- cbind("sex" = rep(sex, nrow(dbsx_ls[[sex]])),
                          "ct" = gsub("\\..*","",rownames(dbsx_ls[[sex]])),
                          dbsx_ls[[sex]])
  }
  dbsx_df <- do.call(rbind, dbsx_ls)
  rownames(dbsx_df) <- NULL
  dbsx_df <- cbind("dbsx"= rep(dbsx, nrow(dbsx_df)), dbsx_df)
  if (dbsx == "DO" | dbsx == "DGN") {
    keep_cols <- c("dbsx", "sex", "ct", "Cluster", "Description", "p.adjust",  "Count")
  } else if (dbsx=="DisGeNET2r_DisGeNET_CURATED") {
    keep_cols <- c("dbsx", "sex", "ct", "Description", "FDR", "Count", "groups")
  } else {
    keep_cols <- c( "dbsx", "sex", "ct", "Term", "Adjusted.P.value", "gene_count", "groups"  )
  }
  dbsx_df <- dbsx_df[, which(colnames(dbsx_df) %in% keep_cols)]
  if (dbsx == "DO" | dbsx == "DGN") {
    dbsx_df <- dbsx_df %>% relocate(Cluster, .after = Count)
  }
  colnames(dbsx_df) <- c("dbsx", "sex", "ct", "term", "adj_pval", "gene_count", "groups")
  return(dbsx_df)
}

# 29. Counts the number of repeated terms
  # Input: main directory where to save the files, the merged dataframe, the order in which plot the groups
  # Return: nothing, saves the files instead

CountDiseases <- function(main_dir, dbsx_all, adj_pval_thresh=0.05, gene_count_thresh=1, which_comp="sex") {
  path <- paste0(main_dir, "/Faceted_Diseases/")
  dir.create(path, recursive = T, showWarnings = F)
  dbsx_all <- subset(dbsx_all, adj_pval <= adj_pval_thresh & gene_count > gene_count_thresh) 
  if (which_comp=="sex") {
    count_sex_df <- list()
    for (sex in unique(dbsx_all$sex)) {
      count_dis <- as.data.frame(table(tolower(dbsx_all[which(dbsx_all$sex==sex), "term"])))
      colnames(count_dis) <- c("term", "freq")
      count_dis <- count_dis[order(count_dis$freq, decreasing = T), ][1:10, ]
      max <- length(unique(dbsx_all[which(dbsx_all$sex==sex), "dbsx"])) * length(unique(dbsx_all[which(dbsx_all$sex==sex), "groups"]))
      count_dis <- cbind(
        "sex" = rep(sex, nrow(count_dis)),
        count_dis,
        "max"= rep(max, nrow(count_dis)))
      count_sex_df <- append(count_sex_df, list(count_dis))
    }
    count_sex_df <- do.call(rbind, count_sex_df)
    count_df <- transform(count_sex_df, perc = freq * 100 / max)
    count_df$term <- as.character(count_df$term)
    write.csv(count_df, paste0(path, "disease_count.csv"))
    #count_df <- complete(count_df, sex, term, fill = list(freq = 0, max = 0, sex_perc =0))
  } else if (which_comp=="sex_ct") {
    count_ct_df <- list()
    for (ct in unique(dbsx_all$ct)) {
      for (sex in unique(dbsx_all$sex)) {
        count_dis <- as.data.frame(table(tolower(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "term"])))
        colnames(count_dis) <- c("term", "freq")
        count_dis <- count_dis[order(count_dis$freq, decreasing = T), ]
        if (nrow(count_dis) > 5) {
          count_dis <- count_dis[1:5,]
        }
        max <- length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "dbsx"])) * length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "groups"]))
        which_dbsx <- vector()
        num_dbsx <- vector()
        for (i in count_dis$term) {
          which_dbsx <- c(which_dbsx, paste(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct & tolower(dbsx_all$term)==i), "dbsx"]), collapse = ", "))
          num_dbsx <- c(num_dbsx, length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct & tolower(dbsx_all$term)==i), "dbsx"])))
        }
        count_dis <- cbind("ct" = rep(ct, nrow(count_dis)),
                           "sex" = rep(sex, nrow(count_dis)),
                           count_dis,
                           "max"= rep(max, nrow(count_dis)),
                           "which_dbsx" = which_dbsx,
                           "num_dbsx" = num_dbsx)
        count_ct_df <- append(count_ct_df, list(count_dis))
      }
    }
    count_ct_df <- do.call(rbind, count_ct_df)
    count_df <- transform(count_ct_df, perc = freq * 100 / max)
    count_df$term <- as.character(count_ct_df$term)
    write.csv(count_df, paste0(path, "disease_count_per_ct.csv"))
  } else if (which_comp=="sex_ct_dbsx") {
    count_ct_df_dbsx <- list()
    for (ct in unique(dbsx_all$ct)) {
      for (sex in unique(dbsx_all$sex)) {
        for (dbsx in unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "dbsx"])) {
          count_dis <- as.data.frame(table(tolower(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct & dbsx_all$dbsx==dbsx), "term"])))
          colnames(count_dis) <- c("term", "freq")
          count_dis <- count_dis[order(count_dis$freq, decreasing = T), ]
          if (nrow(count_dis) > 5) {
            count_dis <- count_dis[1:5,]
          }
          max <- length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct & dbsx_all$dbsx==dbsx), "groups"]))
          count_dis <- cbind("ct" = rep(ct, nrow(count_dis)),
                                  "sex" = rep(sex, nrow(count_dis)),
                                  "which_dbsx" = rep(dbsx, nrow(count_dis)),
                                  count_dis,
                                  "max"= rep(max, nrow(count_dis)))
          count_ct_df_dbsx <- append(count_ct_df_dbsx, list(count_dis))
        }  
      }
    }
    count_ct_df_dbsx <- do.call(rbind, count_ct_df_dbsx)
    count_df <- transform(count_ct_df_dbsx, perc = freq * 100 / max)
    count_df$term <- as.character(count_df$term)
    write.csv(count_df, paste0(path, "disease_count_per_ct_and_dbsx.csv"))
    #count_df <- complete(count_df, ct, sex, which_dbsx, term, fill = list(freq = 0, max = 0, group_perc =0))
  } 
  return(count_df)
}


# 30. Plots the disease enrichemnt results in facets
  # Input: main directory where to save the plots, the merged dataframe, the order in which plot the groups
  # Return: nothing, saves the plots instead

PlotFacetedDB <- function(main_dir, dbsx_all, groups_ordered) {
  dbsx_all <- subset(dbsx_all, adj_pval <= 0.05 & gene_count > 1) 
  dbsx_all$groups <- factor(dbsx_all$groups, groups_ordered[which(groups_ordered %in% unique(dbsx_all$groups))]) 
  dbsx_all <- dbsx_all[order(dbsx_all$groups), ]
  plot_path <- paste0(main_dir, "/Faceted_Diseases/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  for (sex in unique(dbsx_all$sex)) {
    for (ct in unique(dbsx_all[which(dbsx_all$sex==sex), "ct"])) {
      if (nrow(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), ]) > 1) {
        pdf(paste0(plot_path, ct, "_", sex, ".pdf"), width = 25, height = 10)
        print(
          ggplot(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), ], aes(groups, term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(title = paste0(ct, "_", sex), y = "", x = "Groups", size = "Gene count", color = "Adjusted p-value") +
            facet_wrap(~ dbsx, scales = "free", nrow=1) +
            scale_size_continuous(range=c(3, 8)) + 
            scale_y_discrete(labels=function(x) str_wrap(x,width=40)) +
            theme(
              plot.title = element_text(size=14, face="bold", colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=12, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.y = element_blank(),
              axis.text.y = element_text(size=12, colour = "black", vjust = 0.7, hjust=0.5),
              legend.position = "right", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
        )
        dev.off()
      } else {
        print(paste0("The ", ct, " in ", sex, " has no significant terms"))
      }
      
    }
  }
}

# 30. Plots the disease enrichemnt results in facets
  # Input: main directory where to save the plots, the merged dataframe, the order in which plot the groups, which comparison to plot
  # Return: nothing, saves the plots instead

PlotFacetedDBSimplified <- function(main_dir, count_df, which_comp="sex", plot_order="no") {
  plot_path <- paste0(main_dir, "/Faceted_Diseases/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  if (which_comp=="sex") {
    pdf(paste0(plot_path, "top_10_terms_per_sex.pdf"), width = 10, height = 10)
    print(
      ggplot(count_df, aes(reorder(term, -perc), perc, fill=sex)) +
        geom_bar(stat="identity", color="black", position = "dodge") +
        facet_wrap(~sex, scales = "free_x") +
        labs(x="Terms", y="% of groups enriched", fill="Sex") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              strip.text = element_text(size=12, face="bold", colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_blank(),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
      
    )
    dev.off()
  } else if (which_comp=="sex_ct_dbsx") {
    count_df$ct <- factor(count_df$ct, plot_order[which(plot_order %in% count_df$ct)])
    count_df <- count_df[order(count_df$ct),]
    pdf(paste0(plot_path, "top_10_terms_per_sex_ct_dbsx.pdf"), width = 10, height = 15)
    print(
      ggplot(count_df, aes(ct, factor(term, rev(unique(term))),  size = perc , color=which_dbsx)) +
        geom_point() +
        facet_grid(~ sex, scales = "free") +
        labs(x="Cell types", y="Terms", size="% of groups enriched", color="Disease databases") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              strip.text = element_text(size=12, face="bold", colour = "black"),
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
}

# 32. Creates a dataframe with only the gens related to know diseases
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
      deg_presence <- ifelse(ref_df[which(ref_df$group_id==id), "dis_genes"] %in% sex_dfs[which(sex_dfs$id==deg), "gene_id"], "yes", "no")
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

# 33. Plot heatmap with the results of which disease-associated genes are found in the degs
  # Iput: the DEG data frame with the presence of genes-associated genes, the disease group to plot
  # Return: the faceted plot

PlotDisDegGroup <- function(ref_deg, dis_id, groups_ordered) {
  dis_plot <- ggplot(complete(ref_deg[which(ref_deg$disease_group==dis_id),]), aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), dis_gene_id, fill=dis_presence)) +
    geom_tile(color="white") +
    facet_grid(sex ~ ct, scales = "free") +
    labs(x="Groups", y="Disease-associated genes", fill="Genes found", title =dis_id) +
    scale_fill_manual(values = c("yes"="#F8766D",
                                 "no"="#00BFC4"),
                      na.value = "grey",
                      guide = guide_legend(reverse = TRUE)) +
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
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.y = element_blank(),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
}

# 34. Plots the results for the disease-associated genes for all disease groups
  # Input: main directory where to save the plots, the DEG dataframe with the presence of genes-associated genes
  # Return: nothing, saves plots instead

PlotDisDeg <- function(main_dir, ref_deg, ref_df_name, groups_ordered) {
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")  
  dir.create(out_path, showWarnings = F, recursive = T)
  for (dis in unique(ref_deg$disease_group)) {
    print(dis)
    pdf(paste0(out_path, dis, ".pdf"), width = 16)
    print(PlotDisDegGroup(ref_deg, dis, groups_ordered))
    dev.off()
  }
}

# 35. Imports TRANSFAC_and_JASPAR_PWMs results
  # Input: main directory, sub-folders list
  # Return: dataframe containing all significant terms

ImportTJPWMs <- function(main_dir) {
  files_path <- paste0(main_dir, "EnrichR_TRANSFAC_and_JASPAR_PWMs/")
  cts <- list.dirs(files_path, full.names = F, recursive = F)
  tj_ls <- list()
  for (ct_id in cts) {
    tj_files <- list.files(paste0(files_path, ct_id, "/"), pattern = "\\.csv$", full.names = TRUE)
    tj_names <- list.files(paste0(files_path, ct_id, "/"), pattern = "\\.csv$", full.names = F)
    ct_ls <- lapply(tj_files, read.csv, row.names=1)
    names(ct_ls) <- str_remove_all(tj_names, "\\.csv$")
    names(ct_ls) <- str_remove_all(names(ct_ls), "top5_")
    tj_ls <- append(tj_ls, list(ct_ls))
  }
  names(tj_ls) <- cts
  tj_df <- do.call(rbind, unlist(tj_ls, recursive = F))
  tj_groups <- data.frame("groups"=rownames(tj_df))
  rownames(tj_df) <- NULL
  tj_groups <- separate(tj_groups, groups, into=c("ct", "sex", "num"), remove = T, sep="\\.")
  tj_df <- cbind(tj_groups[, c(1:2)], tj_df)
  tj_df <- tj_df[which(tj_df$Adjusted.P.value < 0.05), ]
  tj_df$term_simplified <- str_remove_all(tj_df$Term, " \\(human\\)")
  return(tj_df)
}

# 36. Function to get chromosome number from term
  # Input: the terms as vector
  # Return: the annotated terms with chromosome number

Annot.chr.name <- function(gene.list){
  # define biomart object
  mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "useast")
  Annot_idf <- getBM(attributes = c("hgnc_symbol",
                                    "chromosome_name"),
                     filters = c("hgnc_symbol") ,
                     values=list(gene.list),
                     mart = mart)
  #delete chromosome name with CHR label
  Annot_df <- Annot_idf[!str_detect(Annot_idf$chromosome_name,  "CHR"),]
  return(Annot_df)
}

# 37. Map terms from shared EnrichR_TRANSFAC_and_JASPAR_PWMs against chromosome
  # Input: dataframe with the shared terms, the annotated terms with chromosome number
  # Return: merged dataframe

map_chr <- function(gene_count_filt, Annot_df){
  map_chr_df <- merge(gene_count_filt, Annot_df, by.x= "term", by.y= "hgnc_symbol", all.x=T)
  return(map_chr_df)
}

# 38. Caclulates common EnrichR_TRANSFAC_and_JASPAR_PWMs and saves to CSV with chromosome information
  # Input: main directory where to save the file, the tj dataframe, the percenytage of minimum cts to share each term
  # Return: nothing, saves CSV instead

CalculateSharedTJPWMs <- function(main_dir, tj_df, min_num_ct) {
  out_path <- paste0(main_dir, "EnrichR_TRANSFAC_and_JASPAR_PWMs/")
  sex_tjs_ls <- list()
  thresh <- ceiling(length(unique(tj_df$ct)) * min_num_ct)
  for (sex in c("F", "M")) {
    sex_tjs <- vector()
    for (ct in unique(tj_df$ct)) {
      sex_tjs <- c(sex_tjs, unique(tj_df[which(tj_df$ct==ct & tj_df$sex==sex), "term_simplified"]))
    }
    sex_tjs <- as.data.frame(table(sex_tjs))
    sex_tjs <- cbind("sex" = rep(sex, nrow(sex_tjs)), sex_tjs)
    sex_tjs <- sex_tjs[order(sex_tjs$Freq, decreasing = T), ]
    sex_tjs <- sex_tjs[which(sex_tjs$Freq >= thresh), ]
    sex_tjs_ls <- append(sex_tjs_ls, list(sex_tjs))
  }
  sex_tjs_ls <- do.call(rbind, sex_tjs_ls)
  rownames(sex_tjs_ls) <- NULL
  colnames(sex_tjs_ls) <- c("sex", "term", "ct_count")
  annot_terms <- Annot.chr.name(sex_tjs_ls$term)
  sex_tjs_ls <- map_chr(sex_tjs_ls, annot_terms)
  sex_tjs_ls$chr_simplified <- str_replace_all(sex_tjs_ls$chromosome_name, "\\d+", "Autosome")
  write.csv(sex_tjs_ls, paste0(out_path, "TRANSFAC_and_JASPAR_PWMs_shared_by_", min_num_ct * 100, "%_cts.csv"))
}