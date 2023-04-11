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
library(scales) # to select palette
library(openxlsx) # to import xlsx files
library(RColorBrewer) # for palette

# Sets the EnrichR database for Human genes
setEnrichrSite("Enrichr") 

# Set up the max overlap for label repelling for CNET plots
options(ggrepel.max.overlaps = Inf)

# 3. Compares the GOs of a list of genes
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

# 4. finds the EntrezID for a gene
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

# 5. Compares the KEGG pathways of a list of genes
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
  # Return: enriched DO (formal class compareClusterResult)

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
  # Return: enriched DGN (formal class compareClusterResult)

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




# 14. Compares the DEGs from all groups in a specific ct-sex group
  # Input: main directory where to save the plots, the dataframe containing all DEGs, , the enrichment to be analyzed, 
    # the GO to analyze (if GO is the module), if a minimum threshold for how many genes in each module should be used, 
    # the order in which plot the cts, and if the x-axis labels should be rotated by 90 degrees, the threshold for the adjusted p-value to use
  # Return: nothing, saves plots and CSVs instead

Enrich2ndTrim <- function(main_dir, shared_degs_ls, enrich_module, GO_ont="BP", gene_thresh="no", rotate_x_axis=F, adj_pval_thresh=0.05) {
  out_path <- paste0(main_dir,"Functional_analysis_M_2nd_trim_shared_degs/")
  dir.create(out_path, recursive = T, showWarnings = F)
  if (enrich_module=="GO") {
    try({
      print(paste0("Calculating the ", GO_ont, " results"))
      shared_degs_ls_GO <- compareGO(shared_degs_ls, GO_ont, gene_thresh)
      csv_group_GO <- as.data.frame(shared_degs_ls_GO)
      print("Saving the CSV GO results")
      write.csv(csv_group_GO, paste0(out_path, "GO_", GO_ont, ".csv"))
      shared_degs_ls_GO@compareClusterResult <- shared_degs_ls_GO@compareClusterResult[which(shared_degs_ls_GO@compareClusterResult$p.adjust<=adj_pval_thresh),]
      shared_degs_ls_plot <- dotplot(shared_degs_ls_GO, by = 'count', showCategory=2)
      if (rotate_x_axis) {shared_degs_ls_plot$theme$axis.text.x$angle <- 90}
      png(paste0(out_path, "GO_", GO_ont, ".png"), height = 20, width = 15, units = "cm", res = 300)
      print(shared_degs_ls_plot)
      dev.off()
    })
  } else if (enrich_module=="KEGG") {
    try({
      print("Calculating the KEGG results")
      shared_degs_ls_KEGG <- compareKEGG(shared_degs_ls, gene_thresh)
      csv_group_KEGG <- as.data.frame(shared_degs_ls_KEGG)
      print("Saving the CSV KEGG results")
      write.csv(csv_group_KEGG, paste0(out_path, "KEGG.csv"))
      shared_degs_ls_KEGG@compareClusterResult <- shared_degs_ls_KEGG@compareClusterResult[which(shared_degs_ls_KEGG@compareClusterResult$p.adjust<=adj_pval_thresh),]
      shared_degs_ls_plot <- dotplot(shared_degs_ls_KEGG, by = 'count', showCategory=2)
      if (rotate_x_axis) {shared_degs_ls_plot$theme$axis.text.x$angle <- 90}
      png(paste0(out_path, "KEGG.png"), height = 20, width = 15, units = "cm", res = 300)
      print(shared_degs_ls_plot)
      dev.off()
      png(paste0(out_path, "KEGG_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
      print(cnetplot(setReadable(shared_degs_ls_KEGG, 'org.Hs.eg.db', 'ENTREZID'),cex_label_category = 1.2))
      dev.off()
    })
  } else if (enrich_module=="DO") {
    try({
      print("Calculating the DO results")
      shared_degs_ls_DO <-  compareDO(shared_degs_ls, gene_thresh)
      csv_group_DO <- as.data.frame(shared_degs_ls_DO)
      print("Saving the CSV DO results")
      write.csv(csv_group_DO, paste0(out_path, "DO.csv"))
      shared_degs_ls_DO@compareClusterResult <- shared_degs_ls_DO@compareClusterResult[which(shared_degs_ls_DO@compareClusterResult$p.adjust<=adj_pval_thresh),]
      shared_degs_ls_plot <- dotplot(shared_degs_ls_DO, by = 'count', showCategory=2)
      if (rotate_x_axis) {shared_degs_ls_plot$theme$axis.text.x$angle <- 90}
      png(paste0(out_path, "DO.png"), height = 20, width = 15, units = "cm", res = 300)
      print(shared_degs_ls_plot)
      dev.off()
      png(paste0(out_path, "DO_cnet.png"), height = 15, width = 20, units = "cm", res = 300)
      print(cnetplot(shared_degs_ls_DO, node_label="category",cex_label_category = 1.2))
      dev.off()
    })
  } else if (enrich_module=="DGN") {
    try({
      print("Calculating the DGN results")
      shared_degs_ls_DGN <-  compareDGN(shared_degs_ls, gene_thresh)
      csv_group_DGN <- as.data.frame(shared_degs_ls_DGN)
      print("Saving the CSV DGN results")
      write.csv(csv_group_DGN, paste0(out_path, "DGN.csv"))
      shared_degs_ls_DGN@compareClusterResult <- shared_degs_ls_DGN@compareClusterResult[which(shared_degs_ls_DGN@compareClusterResult$p.adjust<=adj_pval_thresh),]
      shared_degs_ls_plot <- dotplot(shared_degs_ls_DGN, by = 'count', showCategory=2)
      if (rotate_x_axis) {shared_degs_ls_plot$theme$axis.text.x$angle <- 90}
      png(paste0(out_path, "DGN.png"), height = 20, width = 15, units = "cm", res = 300)
      print(shared_degs_ls_plot)
      dev.off()
    })
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

# 21. Calculates the enrichment in the package and database specified for each ct-sex combo across all cts
# Input: main directory where to save the plots, the dataframe containing all DEGs, the package to be used,
# the database, the order in which plot the groups
# Return: nothing, saves plots and CSVs instead

EnrichOtherDBGroup2ndTrim <- function(main_dir, shared_degs_ls, package, dbsx, cts_ordered){
  out_path <- paste0(main_dir,"Functional_analysis_M_2nd_trim_shared_degs/")
  dir.create(out_path, recursive = T, showWarnings = F)
  dbsx_name <- str_replace_all(dbsx, c(" "="_", "\\("="", "\\)"=""))
  filt_flag <- T
  enrich_df <- EnrichId(shared_degs_ls, package, dbsx, "cts")
  top5 <-SelectTop5(enrich_df)
  write.csv(top5, paste0(out_path, dbsx_name, "_top5.csv"))
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
    pdf(paste0(out_path, dbsx_name, ".pdf"))
    print(
      ggplot(filtered_top, aes(factor(cts, cts_ordered[which(cts_ordered %in% cts)]), term, size = gene_count, color = adj_pval)) + 
        geom_point() + 
        guides(size  = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
        scale_color_continuous(low="red", high="blue",guide=guide_colorbar(reverse=T)) +
        labs(title = dbsx_name, y = "", x = "Cell types", size = "Gene count", color = "Adjusted p-value") +
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

