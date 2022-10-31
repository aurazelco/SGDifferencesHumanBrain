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
ct_order <- c(
  "L2_3 EN",               
  "L4 EN",   
  "PLCH1 L4_5 EN", 
  "TSHZ2 L4_5 EN", 
  "L5 EN",       
  "L5_6 EN",       
  "L5b EN",     
  "L6 EN",     
  "pyramidal neuron", 
  "CXCL14 IN",  
  "PVALB IN",                    
  "SST IN",
  "SV2C IN",               
  "VIP IN",  
  "EC", 
  "fibrous astrocyte",
  "protoplasmic astrocyte",
  "OPC", 
  "oligodendrocyte",           
  "microglia"
)

####### 01B_plot_num_genes.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01B_plot_num_genes_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

# NORMAL
num_df_normal <- IntersectDEG(main, sub_disease[3], pval_thresh, FC_thresh, ct_order)
PlotIntDEGs(main, sub_disease[3], num_df_normal[[1]], num_df_normal[[2]], ct_order)

# AD
num_df_AD <- IntersectDEG(main, sub_disease[1], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[1], num_df_AD[[1]], num_df_AD[[2]], ct_order)

# MS
num_df_MS <- IntersectDEG(main, sub_disease[2], pval_thresh, FC_thresh, ct_order)
PlotIntDEGs(main, sub_disease[2], num_df_MS[[1]], num_df_MS[[2]], ct_order)


####### 01C_num_chr.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01C_num_chr_func.R")

# as used in 02A_Fisher
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

# NORMAL
chr_normal <- ProcessCt(main, sub_disease[3])
PlotGeneralHeatmap(main, sub_disease[3], chr_normal, ct_order)
PlotSexHmp(main, sub_disease[3], chr_normal, ct_order)
PlotNumChr(main, sub_disease[3], num_chr_genes, T, ct_order)

# AD
chr_ad <- ProcessCt(main, sub_disease[1])
PlotGeneralHeatmap(main, sub_disease[1], chr_ad, ct_order)
PlotSexHmp(main, sub_disease[1], chr_ad, ct_order)
PlotNumChr(main, sub_disease[1], num_chr_genes, T, ct_order)

# MS
chr_ms <- ProcessCt(main, sub_disease[2])
PlotGeneralHeatmap(main, sub_disease[2], chr_ms, ct_order)
PlotSexHmp(main, sub_disease[2], chr_ms, ct_order)
PlotNumChr(main, sub_disease[2], num_chr_genes, T, ct_order)


####### 01D_Xpar1,2.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01D_Xpar1,2_func.R")

Xpar1 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

# NORMAL
XparCt(main, sub_disease[3], Xpar1_list, Xpar2_list, ct_order)

# AD
XparCt(main, sub_disease[1], Xpar1_list, Xpar2_list, ct_order)

# MS
XparCt(main, sub_disease[2], Xpar1_list, Xpar2_list, ct_order)

####### 01E_CellMarker.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01E_CellMarker_func.R")

ct_list <- c(
  "Astrocyte" ="astrocyte",
  "B cell" = "B",                         
  "Endothelial cell" = "EC",             
  "Glial cell" = "Glia" ,                    
  "Glutamatergic neuron" = "EN",         
  "Interstitial cell" = "IC",              
  "Lake et al.Science.Ex1" = "EN",              
  "Lake et al.Science.Ex2" = "EN",             
  "Lake et al.Science.Ex3" = "EN",              
  "Lake et al.Science.Ex4" = "EN",             
  "Lake et al.Science.Ex5" = "EN",              
  "Lake et al.Science.Ex6" = "EN",             
  "Lake et al.Science.Ex7" = "EN",              
  "Lake et al.Science.Ex8" = "EN",            
  "Lake et al.Science.In1" = "IN",              
  "Lake et al.Science.In2" = "IN",        
  "Lake et al.Science.In3" = "IN",         
  "Lake et al.Science.In4" = "IN",        
  "Lake et al.Science.In5" = "IN",         
  "Lake et al.Science.In6" = "IN",        
  "Lake et al.Science.In7" = "IN",         
  "Lake et al.Science.In8" = "IN",        
  "M1 macrophage" = "Macrophage",                  
  "M2 macrophage" = "Macrophage",                 
  "Macrophage"= "Macrophage",                      
  "Microglial cell" = "Microglia",                
  "Neural progenitor cell" = "NPC",          
  "Neural stem cell" = "NSC",                
  "Neuron" = "Neuron",                          
  "Neutrophil" = "Neutrophil",                    
  "Oligodendrocyte" = "Oligodendrocyte",                
  "Oligodendrocyte precursor cell" = "Oligodendrocyte",
  "Oligodendrocyte progenitor cell" = "Oligodendrocyte",
  "Pericyte"  = "Pericyte",                   
  "Purkinje cell" = "Neuron",                  
  "Stem cell" =    "Stem cell",                 
  "T cell"  = "T",                         
  "T helper2 (Th2) cell" = "T"
)

data_ct <- c("CXCL14 IN" = "IN",
             "EC" = "EC",
             "fibrous astrocyte"  = "astrocyte",
             "L2_3 EN" = "EN", 
             "L4 EN" = "EN",
             "L5 EN" = "EN",
             "L5_6 EN" = "EN",
             "L5b EN" = "EN",
             "L6 EN" = "EN",                
             "microglia" = "Microglia", 
             "oligodendrocyte" =  "Oligodendrocyte",      
             "OPC" = "Oligodendrocyte",                  
             "PLCH1 L4_5 EN" = "EN", 
             "protoplasmic astrocyte" = "astrocyte",
             "PVALB IN"  = "IN",            
             "pyramidal neuron"  = "EN",
             "SST IN" = "IN",   
             "SV2C IN"  = "IN",   
             "TSHZ2 L4_5 EN" = "EN",  
             "VIP IN" = "IN")

extra_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/"

# NORMAL
PlotCMresults(main, sub_disease[3], 
              extra_path, 
              ct_list, data_ct, ct_order)


# AD
PlotCMresults(main, sub_disease[1], 
              extra_path, 
              ct_list, data_ct, ct_order)

# MS
PlotCMresults(main, sub_disease[2], 
              extra_path, 
              ct_list, data_ct, ct_order)



####### 02A_Fisher.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/scripts/DEGs/02A_Fisher_func.R")

# NORMAL
SexChr2(main, sub_disease[3], tot_genes, X_chr_genes, Y_chr_genes)

# AD
SexChr2(main, sub_disease[1], tot_genes, X_chr_genes, Y_chr_genes)

# MS
SexChr2(main, sub_disease[2], tot_genes, X_chr_genes, Y_chr_genes)


####### 02B_ARE_ERE.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02B_ARE_ERE_func.R")

ARE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

# NORMAL
AnalysisARE_ERE(main, sub_disease[3], ARE, EREgene, ct_order)

# AD
AnalysisARE_ERE(main, sub_disease[1], ARE, EREgene, ct_order)

# MS
AnalysisARE_ERE(main, sub_disease[2], ARE, EREgene, ct_order)


####### 02B_ARE_ERE_proj.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02B_ARE_ERE_proj_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

# NORMAL
AnalysisARE_ERE(main, sub_disease[3], pval_thresh, FC_thresh, ARE, EREgene, ct_order)

# AD
AnalysisARE_ERE(main, sub_disease[1], pval_thresh, FC_thresh, ARE, EREgene, ct_order)

# MS
AnalysisARE_ERE(main, sub_disease[2], pval_thresh, FC_thresh, ARE, EREgene, ct_order)

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

SAGD_df <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/SAGD_filt.csv")
SAGD_df[,1] <- NULL

#SAGD <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Sexassociatedgene_Padj0.05_PMID30380119.csv")
#SAGD <- SAGD[, c(2,4)]
#SAGD[which(SAGD[,"Symbol"]==''), "Symbol"]  <- NA
#SAGD <- drop_na(SAGD)
#SAGD$Symbol <- toupper(SAGD$Symbol)
#SAGD$dupl <- paste(SAGD$Symbol, SAGD$Species, sep="_")
#SAGD_drop <- SAGD %>% distinct(dupl, .keep_all = TRUE)
#SAGD_drop <- SAGD_drop[, c(1,2)]

#SAGD_df <- data.frame(unique(SAGD_drop$Symbol))
#colnames(SAGD_df) <- c("gene_name")

#for (sp in SAGD_drop$Species) {
#  SAGD_df[, sp] <- rep(0, length.out = length(SAGD_df$gene))
#}

#for (gene in SAGD_df$gene_name) {
#  species <- SAGD_drop[which(SAGD_drop$Symbol==gene), "Species"]
#  if (length(species) > 0 ) {
#    for (sp in species) {
#      SAGD_df[which(SAGD_df$gene_name==gene), sp] <- 1
#    }
#  }
#}
#names(SAGD_df)[names(SAGD_df) == 'gene'] <- 'gene_name'

#write.csv(SAGD_df, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/SAGD_filt.csv")

# Read all DEGs from healthy patients, regardless of the sex




#library(dplyr)
#norm_deg_counts <- as.data.frame(normal_all_deg %>% 
#                    group_by(cluster) %>%
#                    summarise(tot_degs = length(cluster)))


# all genes commonly expressed in the cts
all_genes <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/tot_genes_ct.csv")
all_genes$X <- NULL
col_factors <- c("disease", "sex", "ct")
all_genes[col_factors] <- lapply(all_genes[col_factors], as.factor) 

# NORMAL
ConservedFractions(main, sub_disease[3], conserved, 4, "Primates", all_genes, ct_order)
ConservedFractions(main, sub_disease[3], SAGD_df, 4, "SAGD",  all_genes, ct_order)

# AD
ConservedFractions(main, sub_disease[1], conserved, 4, "Primates",  all_genes, ct_order)
ConservedFractions(main, sub_disease[1], SAGD_df, 4, "SAGD",  all_genes, ct_order)

# MS
ConservedFractions(main, sub_disease[2], conserved, 4, "Primates",  all_genes, ct_order)
ConservedFractions(main, sub_disease[2], SAGD_df, 4, "SAGD",  all_genes, ct_order)

