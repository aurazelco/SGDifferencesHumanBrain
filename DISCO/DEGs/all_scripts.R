# 01A_generate_DEGs.R
main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)
library(stringr)


disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")
disco_filt@meta.data$proj_sex_disease_ct <- lapply(disco_filt@meta.data$proj_sex_disease_ct, str_replace_all, pattern="/", replacement="_")
Idents(disco_filt) <- "proj_sex_disease_ct"


final_groups <- read.csv("final_filt.csv")
final_groups[, 1] <- NULL

# had to do this because / in ct names were causing issues with dir.create
final_groups$ct <- sapply(final_groups$ct, str_replace_all, pattern="/", replacement="_")
final_groups$idents <- sapply(final_groups$idents, str_replace_all, pattern="/", replacement="_")

col_factors <- c("proj", 
                 "sex",
                 "disease",
                 "ct",
                 "og",
                 "idents"
)

final_groups[col_factors] <- lapply(final_groups[col_factors], as.factor) 


sexes <- levels(final_groups$sex)
for (disease_type in levels(final_groups$disease)) {
  path <- paste0(getwd(), "/", disease_type)
  dir.create(path, showWarnings = FALSE)
  for (ct_type in levels(final_groups$ct)) {
    for (id in levels(final_groups$proj)) {
      id1 <- paste(id, sexes[1], disease_type, ct_type, sep="_")
      id2 <- paste(id, sexes[2], disease_type, ct_type, sep="_")
      if ((id1 %in% levels(final_groups$idents)) & (id2 %in% levels(final_groups$idents))) {
        path <- paste0(getwd(), "/", disease_type, "/", ct_type)
        dir.create(path, showWarnings = FALSE)
        deg1 <- FindMarkers(disco_filt, 
                            ident.1 = id1, 
                            ident.2 = id2,
                            logfc.threshold = 0.25,
                            min.pct = 0.1,
                            only.pos = TRUE)
        write.csv(deg1, paste0(path, "/", id, "_", sexes[1], ".csv"))
        deg2 <- FindMarkers(disco_filt, 
                            ident.1 = id2, 
                            ident.2 = id1,
                            logfc.threshold = 0.25,
                            min.pct = 0.1,
                            only.pos = TRUE)
        write.csv(deg2, paste0(path, "/", id, "_", sexes[2],".csv"))
      }
    }
  }
}

rm(list=ls())

####### GENERAL VARIABLES

main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"
sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

####### 01B_plot_num_genes.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01B_plot_num_genes_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

# NORMAL
num_df_normal <- IntersectDEG(main, sub_disease[3], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[3], num_df_normal[[1]], num_df_normal[[2]])

# AD
num_df_AD <- IntersectDEG(main, sub_disease[1], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[1], num_df_AD[[1]], num_df_AD[[2]])

####### 01C_num_chr.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01C_num_chr_func.R")

# as used in 02A_Fisher
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

# NORMAL
chr_normal <- ProcessCt(main, sub_disease[3])
PlotSexHmp(main, sub_disease[3], chr_normal)
PlotNumChr(main, sub_disease[3], num_chr_genes, T)

#AD
chr_ad <- ProcessCt(main, sub_disease[1])
PlotSexHmp(main, sub_disease[1], chr_ad)
PlotNumChr(main, sub_disease[1], num_chr_genes, T)

####### 01D_Xpar1,2.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01D_Xpar1,2_func.R")

Xpar1 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

# NORMAL
normal_df <- XparCt(main, sub_disease[3], Xpar1_list, Xpar2_list)

#AD
ad_df <- XparCt(main, sub_disease[1], Xpar1_list, Xpar2_list)

####### 02A_Fisher.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02A_Fisher_func.R")

# NORMAL
SexChr2(main, sub_disease[3], tot_genes, X_chr_genes, Y_chr_genes)

#AD
SexChr2(main, sub_disease[1], tot_genes, X_chr_genes, Y_chr_genes)

####### 02B_ARE_ERE.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02B_ARE_ERE_func.R")

ARE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

# NORMAL
AnalysisARE_ERE(main, sub_disease[3], ARE, EREgene)

#AD
AnalysisARE_ERE(main, sub_disease[1], ARE, EREgene)

####### 02C_Conservation.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02C_Conservation_func.R")

# CONSERVATION ACROSS PRIMATES
conserved <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/mart_export.txt",
                      sep = '\t',
                      header = TRUE,
                      fill = TRUE)

colnames(conserved) <- c("gene_id",
                         "gene_name",
                         "source",
                         "gene_GC_cont",
                         "Bolivian_Squirrel_Monkey",
                         "Chimpanzee",             
                         "Gorilla",                 
                         "Gibbon",                 
                         "Olive_Baboon",            
                         "Macaque")

for (sp in seq(5,10)) {
  conserved[which(startsWith(conserved[,sp],"EN")),sp] <- 1
  conserved[which(conserved[,sp]==''), sp] <- 0
  conserved[,sp] <- as.numeric(conserved[,sp])
}
conserved <- conserved %>% distinct(gene_name, .keep_all = TRUE)

# SAGD CONSERVATION
SAGD <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Sexassociatedgene_Padj0.05_PMID30380119.csv")
SAGD <- SAGD[, c(2,4)]
SAGD[which(SAGD[,"Symbol"]==''), "Symbol"]  <- NA
SAGD <- drop_na(SAGD)
SAGD$Symbol <- toupper(SAGD$Symbol)
SAGD$dupl <- paste(SAGD$Symbol, SAGD$Species, sep="_")
SAGD_drop <- SAGD %>% distinct(dupl, .keep_all = TRUE)
SAGD_drop <- SAGD_drop[, c(1,2)]

SAGD_df <- data.frame(unique(SAGD_drop$Symbol))
colnames(SAGD_df) <- c("gene_name")

for (sp in SAGD_drop$Species) {
  SAGD_df[, sp] <- rep(0, length.out = length(SAGD_df$gene))
}

for (gene in SAGD_df$gene_name) {
  species <- SAGD_drop[which(SAGD_drop$Symbol==gene), "Species"]
  if (length(species) > 0 ) {
    for (sp in species) {
      SAGD_df[which(SAGD_df$gene_name==gene), sp] <- 1
    }
  }
}

# NORMAL
ConservedFractions(main, sub_disease[3], conserved, 4, "conserved")
ConservedFractions(main, sub_disease[3], SAGD_df, 4, "SAGD")


#AD
ConservedFractions(main, sub_disease[1], conserved, 4, "conserved")
ConservedFractions(main, sub_disease[1], SAGD_df, 4, "SAGD")
