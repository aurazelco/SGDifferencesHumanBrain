# Author: Aura Zelco
# Brief description:
  # This script is used to run the enrichment analysis on multiple databases/terms, all using the DEGs from the projects
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all conditions
  # 3. For each enrichments, compares F v M in each ct-condition combo and in each sex-ct combo across conditions
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

# 2. Function to filter out ns genes and too low FC, and order based on FC 
  # Input: dataframe of DEGs
  # Return: gene list of significant genes as data.frame

Filter_gene <- function( order.gene.df, pval, FC) {
  logFC <- log2(1/FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  
  #If there are sig genes add index number for each sig gene
  return(data.frame("Genes"=rownames(gene.sig)))
}


# 3. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
  # Input: directory where to find ct sub-folders, if UCSC or not, list of projects ids, the individual project id to look for, 
    # the threshold for p-value and Fc if unfiltered data are imported, file extension, where to find row-names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCt <- function(main_dir, UCSC_flag="no", individual_projs=vector(), single_proj="", pval, FC, ext, row_col) {
  if (length(individual_projs)==0) {
    path <- paste0(main_dir, "/01B_num_DEGs")
  } else {
    path <- paste0(main_dir, "/01A_DEGs")
  }
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportFiles(paste(path, sub_ct[ct], sep="/"))
    if (length(individual_projs)>0) {
      if (length(grep(single_proj, names(deg))) > 0) {
        deg <- deg[grepl(single_proj, names(deg),fixed = T)]
      } else {
        next
      }
    }
    for (i in names(deg)) {
      if (UCSC_flag=="no") {
        if (length(individual_projs)==0) {
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
          deg_ct <- Filter_gene(deg[[i]], pval, FC)
          rownames(deg_ct) <- NULL
          if (grepl("F", i, fixed=TRUE)){
            df_F <- append(df_F, list(deg_ct))
            names_F <- c(names_F, sub_ct[ct])
          } else {
            df_M <- append(df_M, list(deg_ct))
            names_M <- c(names_M, sub_ct[ct])
          }
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
  df_F <- df_F[lapply(df_F,length)>0]
  df_M <- df_M[lapply(df_M,length)>0]
  if (length(df_F) != 0 & length(df_F) != 0) {
    return(list("F"=df_F, "M"=df_M))
  } else {
    return("empty")
  }
}

# 4. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, sub-folders list, if UCSC or not, list of projects ids, , the threshold for p-value and Fc if unfiltered data are imported
  # Return: list of condition lists, each containing ct lists divided in F and M

ImportDataset <- function(main_dir, folder_list, UCSC_flag="no", individual_projs=vector(), pval, FC) {
  ds_list <- list()
  ct_list <- vector()
  group_names <- vector()
  for (folder in folder_list) {
    if (UCSC_flag=="no") {
      if (length(individual_projs)==0) {
        ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder), pval, FC)))
        ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
        group_names <- c(group_names, folder)
      } else {
        for (single_proj in individual_projs) {
          single_proj_list <- list(ImportCt(paste0(main_dir, folder), individual_projs = individual_projs, single_proj = single_proj, pval, FC))
          if (single_proj_list!="empty") {
            ds_list <- append(ds_list, single_proj_list)
            ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/01A_DEGs"), recursive=FALSE, full.names = FALSE))
            group_names <- c(group_names, paste(folder, single_proj, sep = "_"))
          }
        }
      }
    } else {
      if (length(individual_projs)==0) {
        ds_list <- append(ds_list, list(ImportCt(paste0(main_dir, folder, "/outputs"), UCSC_flag)))
        ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/outputs/01B_num_DEGs"), recursive=FALSE, full.names = FALSE))
        group_names <- c(group_names, folder)
      } else {
        for (single_proj in individual_projs) {
          single_proj_list <- list(ImportCt(paste0(main_dir, folder), UCSC_flag, individual_projs = individual_projs, single_proj = single_proj, pval, FC))
          if (single_proj_list!="empty") {
            ds_list <- append(ds_list, single_proj_list)
            ct_list <-c(ct_list, list.dirs(paste0(main_dir, folder, "/outputs/01A_DEGs"), recursive=FALSE, full.names = FALSE))
            group_names <- c(group_names, paste(folder, single_proj, sep = "_"))
          }
          
        }
      }
    }
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(list("genes"=ds_list, "ct"=unique(ct_list)))
}

# 5. Groups cts according to common annotation
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
        sex_df <- separate(sex_df, groups, into=c("condition", "sex", "ct", "gene_num"), sep="\\.")
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

# 6. Creates a list, each element the gene list from each condition present for a specific ct and sex
  # Input: the dataframe, the ct and sex to be analyzed, the order in which the groups should be plotted
  # Return: the list of gene lists, from a specific ct and sex, for each condition

ExtractSexCt <- function(sex_df, ct, sex, condition_ordered) {
  sex_ct <- sex_df[which(sex_df$common_annot==ct & sex_df$sex==sex), ]
  sex_ct <- split(sex_ct$gene_id, sex_ct$condition)
  sex_order <- condition_ordered[which(condition_ordered %in% names(sex_ct))]
  sex_ct <- sex_ct[sex_order]
  return(sex_ct)
}

# 7. Compares the GOs of a list of genes
  # Input: list of genes to be compared, which GO (BP, MO or CC),  
    # if a minimum threshold for how many genes in each module should be used
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
  return(gene.df$ENTREZID)})
  length(sex_list_kegg)
  names(sex_list_kegg) <- names(sex_list)
  print(sex_list_kegg)
  enrich <- compareCluster(geneCluster  = sex_list_kegg, 
                           fun          = "enrichKEGG",
                           pAdjustMethod= "BH",
                           pvalueCutoff = 0.05,
                           qvalueCutoff = 0.05)
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
  for (cond in unique(sex_df$condition)) {
    print(cond)
    cond_path <- paste0(out_path, cond, "/")
    dir.create(cond_path, recursive = T, showWarnings = F)
    cond_df <- sex_df[which(sex_df$condition==cond), ]
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


# 13. Compares the DEGs from all condition in a specific ct-sex group
  # Input: main directory where to save the plots, the dataframe containing all DEGs, , the enrichment to be analyzed, 
    # the GO to analyze (if GO is the module), if a minimum threshold for how many genes in each module should be used, 
    # the order in which plot the conditions, and if the x-axis labels should be rotated by 90 degrees
  # Return: nothing, saves plots and CSVs instead

EnrichCondition <- function(main_dir, sex_df, enrich_module, GO_ont="BP", gene_thresh="no", condition_ordered, rotate_x_axis=F) {
  out_path <- paste0(main_comparison, enrich_module, "_comparison_cts/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (ct in unique(sex_df$common_annot)) {
    print(ct)
    ct_path <- paste0(out_path, ct, "/")
    dir.create(ct_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      sex_ct <- ExtractSexCt(sex_df, ct, sex, condition_ordered)
      if (enrich_module=="GO") {
        try({
          print(paste0("Calculating the ", GO_ont, " results for ", sex))
          sex_cond_GO <-  compareGO(sex_ct, GO_ont, gene_thresh)
          sex_cond_plot <- dotplot(sex_cond_GO, by = 'count', title=sex)
          if (rotate_x_axis) {sex_cond_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "GO_", GO_ont, "_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_cond_plot)
          dev.off()
          csv_cond_GO <- as.data.frame(sex_cond_GO)
          print(paste0("Saving the CSV GO results for ", sex))
          write.csv(csv_cond_GO, paste0(ct_path, "GO_", GO_ont, "_", sex, ".csv"))
        })
      } else if (enrich_module=="KEGG") {
        try({
          print(paste0("Calculating the KEGG results for ", sex))
          sex_cond_KEGG <-  compareKEGG(sex_ct, gene_thresh)
          sex_cond_plot <- dotplot(sex_cond_KEGG, by = 'count', title=sex)
          if (rotate_x_axis) {sex_cond_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "KEGG_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_cond_plot)
          dev.off()
          png(paste0(ct_path, "KEGG_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(setReadable(sex_cond_KEGG, 'org.Hs.eg.db', 'ENTREZID'),cex_label_category = 1.2))
          dev.off()
          csv_cond_KEGG <- as.data.frame(sex_cond_KEGG)
          print(paste0("Saving the CSV KEGG results for ", sex))
          write.csv(csv_cond_KEGG, paste0(ct_path, "KEGG_", sex, ".csv"))
        })
      } else if (enrich_module=="DO") {
        try({
          print(paste0("Calculating the DO results for ", sex))
          sex_cond_DO <-  compareDO(sex_ct, gene_thresh)
          sex_cond_plot <- dotplot(sex_cond_DO, by = 'count', title=sex)
          if (rotate_x_axis) {sex_cond_plot$theme$axis.text.x$angle <- 90}
          png(paste0(ct_path, "DO_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
          print(sex_cond_plot)
          dev.off()
          png(paste0(ct_path, "DO_", sex, "_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
          print(cnetplot(sex_cond_DO, node_label="category",cex_label_category = 1.2))
          dev.off()
          csv_cond_DO <- as.data.frame(sex_cond_DO)
          print(paste0("Saving the CSV DO results for ", sex))
          write.csv(csv_cond_DO, paste0(ct_path, "DO_", sex, ".csv"))
        })
      } else if (enrich_module=="DGN") {
        try({
        print(paste0("Calculating the DGN results for ", sex))
        sex_cond_DGN <-  compareDGN(sex_ct, gene_thresh)
        sex_cond_plot <- dotplot(sex_cond_DGN, by = 'count', title=sex)
        if (rotate_x_axis) {sex_cond_plot$theme$axis.text.x$angle <- 90}
        png(paste0(ct_path, "DGN_", sex, ".png"), height = 20, width = 15, units = "cm", res = 300)
        print(sex_cond_plot)
        dev.off()
        csv_cond_DGN <- as.data.frame(sex_cond_DGN)
        print(paste0("Saving the CSV DGN results for ", sex))
        write.csv(csv_cond_DGN, paste0(ct_path, "DGN_", sex, ".csv"))
        })
      }
    }
  }
}


# 14. Searches the selected database in enrichR
  # Input: the gene list, the database to look into
  # Return: dataframe with the retrieved information

EnrichR_fun <- function(gene_ls, dbsx){
  enrichR_obj <- enrichr(gene_ls, dbsx)
  enrichR_df <- enrichR_obj[[1]]
  gene_count <-  unlist(lapply(enrichR_df$Overlap, function(x) { unlist(strsplit(x, "/"))[1] } ))
  enrichR_df$gene_count <- gene_count
  return(enrichR_df)
}

# 15. Searches in the DisGeNET database using disgenet2r
  # Input: the gene list
  # Return: dataframe with the retrieved information

EnrichDisgenet2r_fun <- function(gene_ls){
  dgn_res <-disease_enrichment( entities = gene_ls, 
                                vocabulary = "HGNC", 
                                database = "CURATED")
  dgn_res <- dgn_res@qresult
  return(dgn_res)
}

# 16. Searches in the selected packages and databases for disease-enrichment
  # Input:  list of genes to be compared, the dtabase to look into, the package to use
  # Return: list of enriched terms for each condition

EnrichCt <- function(sex_ct, package, dbsx){
  if(package == 'EnrichR'){
    enrich_ls <- lapply(sex_ct, function(x) EnrichR_fun(x, dbsx))
  }
  else if(package == 'DisGeNET2r'){
    enrich_ls <- lapply(sex_ct, function(x) EnrichDisgenet2r_fun(x))
  }
  for (cond in names(enrich_ls)) {
    cond_df <- enrich_ls[[cond]]
    if  (!is.null(dim(cond_df)) ){
      cond_df$condition <- rep(cond, nrow(cond_df))
      enrich_ls[[cond]] <- cond_df
    }
  }
  return(enrich_ls)
}

# 17. Selects the top 5 enriched terms to combine in one df for plotting
  # Input: the enriched list
  # Return: the filtered list of terms

SelectTop5 <- function(enrich_df){
  head_df <- lapply(enrich_df, function(x) head(x, 5))
  combind_df <- Reduce(rbind, head_df)
  # Add gene ratio column
  #combind_df$gene_ratio <- sapply(combind_df$Overlap, function(x) eval(parse(text = x)) )
  return(combind_df)
}

# 18. Calculates the enrichment in the package and database specified for each ct-sex combo across all cts
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
    # the database, the order in which plot the conditions
  # Return: nothing, saves plots and CSVs instead

EnrichOtherDB <- function(main_dir, sex_df, package, dbsx, condition_ordered){
  dbsx_path <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  out_path <- paste0(main_dir, package, "_", dbsx_path, "/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (ct in unique(sex_df$common_annot)) {
    print(ct)
    ct_path <- paste0(out_path, ct, "/")
    dir.create(ct_path, recursive = T, showWarnings = F)
    for (sex in c("F", "M")) {
      filt_flag <- T
      sex_ct <- ExtractSexCt(sex_df, ct, sex, condition_ordered)
      enrich_df <- EnrichCt(sex_ct, package, dbsx)
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
        filtered_top <- filtered_top[, c("condition", y_var, size_var, color_var)]
        colnames(filtered_top) <- c("condition", "term", "gene_count", "adj_pval")
        pdf(paste0(ct_path, sex, ".pdf"))
        print(
          ggplot(filtered_top, aes(condition, term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(title = sex, y = "", x = "Developmental Conditions", size = "Gene count", color = "Adjusted p-value") +
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

# 19. Calculates the enrichment in the package and database specified comparing for each ct-condition combo the sexes
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
    # the database, the order in which plot the conditions
  # Return: nothing, saves plots and CSVs instead

EnrichOtherDBFvM <- function(main_dir, sex_df, package, dbsx, condition_ordered){
  dbsx_path <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  out_path <- paste0(main_dir, package, "_", dbsx_path, "_ct_sex/")
  dir.create(out_path, recursive = T, showWarnings = F)
  for (cond in unique(sex_df$condition)) {
    print(cond)
    cond_df <- sex_df[which(sex_df$condition==cond), ]
    for (ct in unique(cond_df$common_annot)) {
      cond_path <- paste0(out_path, cond, "/", ct, "/")
      dir.create(cond_path, recursive = T, showWarnings = F)
      filt_flag <- T
      print(ct)
      f_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="F"), "gene_id"]
      m_genes <- cond_df[which(cond_df$common_annot==ct & cond_df$sex=="M"), "gene_id"]
      sex_comp <- list("F"=f_genes, "M"=m_genes)
      enrich_df <- EnrichCt(sex_comp, package, dbsx)
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
        filtered_top <- filtered_top[, c("condition", y_var, size_var, color_var)]
        colnames(filtered_top) <- c("condition", "term", "gene_count", "adj_pval")
        pdf(paste0(cond_path, "sexes.pdf"))
        print(
          ggplot(filtered_top, aes(condition, term, size = gene_count, color = adj_pval)) + 
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

# 20. Plot the presence heatmap
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotHmpRef <- function(ref_presence_df, ref_ct_id, plot_titles) {
  ref_plot <- ggplot(ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ], aes(cond_ct, gene_ids, fill=presence)) +
    geom_tile() +
    facet_grid(condition ~ sex, scales = "free") +
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

# 21. Plot the presence number of genes
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotBarPlotRef <- function(ref_presence_df, ref_ct_id, plot_titles) {
  ref_plot <- ggplot(ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ], aes(cond_ct, fill = presence)) +
    geom_bar(position = "stack") +
    facet_grid(sex ~ condition, scales = "free") +
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

# 22. Calculates the % of known markers in the DEGs
  # Input: the presence df, the ct to plot, the gene lists
  # Return: the percent df

RefPerc <- function(ref_presence_df, ref_ct_id, sex_df) {
  pos_markers <- ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ]
  tot_genes <- vector()
  tot_names <- vector()
  num_pos <- vector()
  for (id in unique(pos_markers$condition)) {
      for (ct in unique(pos_markers[which(pos_markers$condition==id), "cond_ct"])) {
        tot_names <- c(tot_names, paste(id, ct, "F", sep = "/"), paste(id, ct, "M", sep = "/"))
        num_pos <- c(num_pos, length(pos_markers[which(pos_markers$condition==id & pos_markers$cond_ct==ct & pos_markers$sex=="F" & pos_markers$presence=="Yes"), "gene_ids"]))
        num_pos <- c(num_pos, length(pos_markers[which(pos_markers$condition==id & pos_markers$cond_ct==ct & pos_markers$sex=="M" & pos_markers$presence=="Yes"), "gene_ids"]))
        tot_genes <- c(tot_genes, length(sex_df[which(sex_df$condition==id & sex_df$common_annot==ct & sex_df$sex=="F"), "gene_id"]))
        tot_genes <- c(tot_genes, length(sex_df[which(sex_df$condition==id & sex_df$common_annot==ct & sex_df$sex=="M"), "gene_id"]))
      }
    }
  ref_perc <- data.frame(tot_names, num_pos, tot_genes)
  ref_perc <- separate(ref_perc, tot_names, into = c("condition", "ct", "sex"), sep = "/")
  ref_perc$perc <- ref_perc$num_pos * 100 / ref_perc$tot_genes
  return(ref_perc)
}

# 23. Plot the presence number of genes as percentage
  # Input: the presence df, the ct to plot
  # Return: the plot

PlotBarPlotRefPerc <- function(ref_perc, ref_ct_id, plot_titles) {
  ref_plot <- ggplot(ref_perc, aes(ct, perc, fill=condition)) +
    geom_bar(stat="identity", color="black") +
    facet_grid(sex ~ condition, scales = "free") +
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


# 24. Plots if gebes from a reference df are found or not in the DEGs
  # Input: main directory where to save the plots, the dataframe containing all DEGs, the reference df, 
    # the order in which plot the conditions (and which conditions to plot), the vector to use for plot titles
  # Return: nothing, saves plot instead

PlotRefCt <- function(main_dir, sex_df, ref_df, condition_ordered, ref_df_name="ref", plot_titles){
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")
  dir.create(out_path, recursive = T, showWarnings = F)
  sex_df_filt <- sex_df[which(sex_df$condition %in% condition_ordered), ]
  presence <- vector()
  ids <- vector()
  gene_ids <- vector()
  for (sex_id in unique(sex_df_filt$sex)) {
    for (ct in unique(ref_df$Celltype)) {
      ref_genes <- ref_df[which(ref_df$Celltype==ct), "gene"]
      for (cond in unique(sex_df_filt[which(sex_df_filt$sex==sex_id), "condition"])) {
        for (ct_id in unique(sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$condition==cond), "common_annot"])) {
          presence <- c(presence, 
                        ifelse(ref_df[which(ref_df$Celltype==ct), "gene"] %in% sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$condition==cond & sex_df_filt$common_annot==ct_id), "gene_id"],
                               "Yes", "No"))
          ids <- c(ids, 
                   rep(paste(sex_id, ct, cond, ct_id, sep = "/"), length(ref_genes)))
          gene_ids <- c(gene_ids, ref_genes)
          
        }
      }
    }
  }
  ref_presence_df <- data.frame(ids, gene_ids, presence)
  ref_presence_df <- separate(ref_presence_df, ids, into=c("sex", "ref_ct", "condition", "cond_ct"), sep = "/")
  ref_presence_df$condition <- factor(ref_presence_df$condition, condition_order)
  ref_presence_df <- ref_presence_df[order(ref_presence_df$condition), ]
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
    print(PlotBarPlotRefPerc(ref_perc, ref_ct_id, plot_titles))
    dev.off()
  }
}


# 25. Import results from disease-related enrichments
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
    keep_cols <- c("dbsx", "sex", "ct", "Description", "FDR", "Count", "condition")
  } else {
    keep_cols <- c( "dbsx", "sex", "ct", "Term", "Adjusted.P.value", "gene_count","condition"  )
  }
  dbsx_df <- dbsx_df[, which(colnames(dbsx_df) %in% keep_cols)]
  if (dbsx == "DO" | dbsx == "DGN") {
    dbsx_df <- dbsx_df %>% relocate(Cluster, .after = Count)
  }
  colnames(dbsx_df) <- c("dbsx", "sex", "ct", "term", "adj_pval", "gene_count", "condition")
  return(dbsx_df)
}

# 26. Counts the number of repeated terms
  # Input: main directory where to save the files, the merged dataframe, the order in which plot the conditions
  # Return: nothing, saves the files instead

CountDiseases <- function(main_dir, dbsx_all) {
  path <- paste0(main_dir, "/Faceted_Diseases/")
  dir.create(path, recursive = T, showWarnings = F)
  dbsx_all <- subset(dbsx_all, adj_pval <= 0.05 & gene_count > 1) 
  count_ct_df <- list()
  for (ct in unique(dbsx_all$ct)) {
    for (sex in unique(dbsx_all$sex)) {
      count_dis <- as.data.frame(table(tolower(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "term"])))
      count_dis <- count_dis[order(count_dis$Freq, decreasing = T), ][1:10, ]
      max <- length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "dbsx"])) * length(unique(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), "condition"]))
      count_dis <- cbind("ct" = rep(ct, nrow(count_dis)),
                         "sex" = rep(sex, nrow(count_dis)),
                         count_dis,
                         "max"= rep(max, nrow(count_dis)))
      count_ct_df <- append(count_ct_df, list(count_dis))
    }
  }
  count_ct_df <- do.call(rbind, count_ct_df)
  write.csv(count_ct_df, paste0(path, "disease_count_per_ct.csv"))
  count_sex_df <- list()
  for (sex in unique(dbsx_all$sex)) {
    count_dis <- as.data.frame(table(tolower(dbsx_all[which(dbsx_all$sex==sex), "term"])))
    count_dis <- count_dis[order(count_dis$Freq, decreasing = T), ][1:20, ]
    max <- length(unique(dbsx_all[which(dbsx_all$sex==sex), "dbsx"])) * length(unique(dbsx_all[which(dbsx_all$sex==sex), "condition"]))
    count_dis <- cbind(
                       "sex" = rep(sex, nrow(count_dis)),
                       count_dis,
                       "max"= rep(max, nrow(count_dis)))
    count_sex_df <- append(count_sex_df, list(count_dis))
  }
  count_sex_df <- do.call(rbind, count_sex_df)
  write.csv(count_sex_df, paste0(path, "disease_count.csv"))
}


# 27. Plots the disease enrichemnt results in facets
  # Input: main directory where to save the plots, the merged dataframe, the order in which plot the conditions
  # Return: nothing, saves the plots instead

PlotFacetedDB <- function(main_dir, dbsx_all, condition_ordered) {
  dbsx_all <- subset(dbsx_all, adj_pval <= 0.05 & gene_count > 1) 
  dbsx_all$condition <- factor(dbsx_all$condition, condition_ordered[which(condition_ordered %in% unique(dbsx_all$condition))]) 
  dbsx_all <- dbsx_all[order(dbsx_all$condition), ]
  plot_path <- paste0(main_dir, "/Faceted_Diseases/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  for (sex in unique(dbsx_all$sex)) {
    for (ct in unique(dbsx_all[which(dbsx_all$sex==sex), "ct"])) {
      if (nrow(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), ]) > 1) {
        pdf(paste0(plot_path, ct, "_", sex, ".pdf"), width = 25, height = 10)
        print(
          ggplot(dbsx_all[which(dbsx_all$sex==sex & dbsx_all$ct==ct), ], aes(condition, term, size = gene_count, color = adj_pval)) + 
            geom_point() + 
            guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
            scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
            labs(title = paste0(ct, "_", sex), y = "", x = "Developmental Conditions", size = "Gene count", color = "Adjusted p-value") +
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

# 28. Creates a dataframe with only the gens related to know diseases
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
  sex_dfs$id <- paste(sex_dfs$condition, sex_dfs$sex, sex_dfs$common_annot, sep = "/")
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
  ref_deg <- separate(ref_deg, deg_ids, into = c("condition", "sex", "ct"), sep = "/")
  ref_deg$dis_gene_id <- paste(ref_deg$disease, ref_deg$genes_ids, sep =  " - ")
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")  
  dir.create(out_path, showWarnings = F, recursive = T)
  write.csv(ref_deg, paste0(out_path, "disease_genes_in_degs.csv"))
  return(ref_deg)
}

# 29. Plot heatmap with the results of which disease-associated genes are found in the degs
  # Iput: the DEG data frame with the presence of genes-associated genes, the disease group to plot
  # Return: the faceted plot

PlotDisDegGroup <- function(ref_deg, dis_id) {
  dis_plot <- ggplot(complete(ref_deg[which(ref_deg$disease_group==dis_id),]), aes(factor(condition, condition_order), dis_gene_id, fill=dis_presence)) +
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

# 30. Plots the results for the disease-associated genes for all disease groups
  # Input: main directory where to save the plots, the DEG dataframe with the presence of genes-associated genes
  # Return: nothing, saves plots instead

PlotDisDeg <- function(main_dir, ref_deg, ref_df_name) {
  out_path <- paste0(main_dir, "Hmp_", ref_df_name, "/")  
  dir.create(out_path, showWarnings = F, recursive = T)
  for (dis in unique(ref_deg$disease_group)) {
    print(dis)
    pdf(paste0(out_path, dis, ".pdf"), width = 16)
    print(PlotDisDegGroup(ref_deg, dis))
    dev.off()
  }
}