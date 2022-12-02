# 01A_generate_DEGs.R

library(Seurat)
library(stringr)

main <- "/Home/ii/auraz/data/UCSC/outputs/DEGs/"

dir.create(main, recursive = T, showWarnings = F)

trim_2nd <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_2nd_trimester.rds")
trim_2nd@meta.data$sex_ct <- paste(trim_2nd@meta.data$sex, trim_2nd@meta.data$cluster_final, sep="_")
Idents(trim_2nd) <- "sex_ct"

min_cells <- 100

main_deg <- paste0(main, "Velmeshev_2022_2nd_trimester/outputs/")
dir.create(main_deg, recursive = T, showWarnings = F)

final_groups <- read.csv(paste0(main_deg, "final_filt_", min_cells, ".csv"))
sexes <- unique(final_groups$sex)
sexes <- c("F"="Female", "M"="Male")


path <- paste0(main_deg, "01A_DEGs/")
dir.create(path, showWarnings = FALSE, recursive = T)
for (ct_type in unique(final_groups$ct)) {
  id1 <- paste(sexes[[1]], ct_type, sep="_")
  id2 <- paste(sexes[[2]], ct_type, sep="_")
  if ((id1 %in% unique(final_groups$idents)) & (id2 %in% unique(final_groups$idents))) {
    path_ct <- paste0(path, "/", ct_type, "/")
    dir.create(path_ct, showWarnings = FALSE)
    deg1 <- FindMarkers(trim_2nd, 
                              ident.1 = id1, 
                              ident.2 = id2,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
    write.csv(deg1, paste0(path_ct, names(sexes[1]), ".csv"))
    deg2 <- FindMarkers(trim_2nd, 
                              ident.1 = id2, 
                              ident.2 = id1,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
    write.csv(deg2, paste0(path_ct, names(sexes[2]),".csv"))
  }
}


#write.csv(final_groups, paste0(main, "final_filt_", min_cells, ".csv"),row.names = F)

rm(list=ls())

####### GENERAL VARIABLES

main <- "/Home/ii/auraz/data/UCSC/outputs/DEGs/Velmeshev_2022_2nd_trimester/outputs/"
main_local <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2nd_trimester/outputs/"

ct_order <- c(
  "Dorsal progenitors",
  "Ventral progenitors",
  "Excitatory neurons",    
  "Interneurons",                      
  "OPCs",  
  "Astrocytes",
  "Microglia",   
  "Vascular cells",  
  "Unknown"
)

####### 01B_plot_num_genes.R

source("/Home/ii/auraz/scripts/DEGs/01B_plot_num_genes_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

# NORMAL
CountDEG(main, pval_thresh, FC_thresh, ct_order)

####### 01C_num_chr.R -> not on the server

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/DEGs/01C_num_chr_func.R")

# as used in 02A_Fisher
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

chr_2nd_trim <- ProcessCt(main_local)
PlotGeneralHeatmap(main_local, chr_2nd_trim, ct_order)
# doe snot work -> need to check later why not, but these plots are non-essentials
#PlotSexHmp(main_local, chr_2nd_trim, ct_order)


####### 01D_Xpar1,2.R

source("/Home/ii/auraz/scripts/DEGs/01D_Xpar1,2_func.R")

Xpar1 <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

XparCt(main, Xpar1_list, Xpar2_list, ct_order)

####### 02A_Fisher.R  -> not on the server

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/DEGs/02A_Fisher_func.R")

SexChr2(main_local, tot_genes, X_chr_genes, Y_chr_genes)
PlotNumChr(main_local, num_chr_genes, ct_order, T)

####### 02B_ARE_ERE.R

source("/Home/ii/auraz/scripts/DEGs/02B_ARE_ERE_func.R")

ARE <- read_excel("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

AnalysisARE_ERE(main, ARE, EREgene, ct_order)

####### 02C_Conservation.R

source("/Home/ii/auraz/scripts/DEGs/02C_Conservation_func.R")

# CONSERVATION ACROSS PRIMATES
conserved <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/mart_export.txt",
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

SAGD_df <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/SAGD_filt.csv")
SAGD_df[,1] <- NULL

#SAGD <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/Sexassociatedgene_Padj0.05_PMID30380119.csv")
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

#write.csv(SAGD_df, "/Home/ii/auraz/data/UCSC/outputs/DEGs/extra_files/SAGD_filt.csv")

# Read all DEGs from healthy patients, regardless of the sex




#library(dplyr)
#norm_deg_counts <- as.data.frame(normal_all_deg %>% 
#                    group_by(cluster) %>%
#                    summarise(tot_degs = length(cluster)))


# all genes commonly expressed in the cts
all_genes <- read.csv("/Home/ii/auraz/data/UCSC/outputs/DEGs/Velmeshev_2022_2nd_trimester/tot_genes_ct_Velmeshev_2022_2nd_trimester.csv")
all_genes$X <- NULL
col_factors <- c("sex", "ct")
all_genes[col_factors] <- lapply(all_genes[col_factors], as.factor) 


ConservedFractions(main, conserved, 4, "Primates", all_genes, ct_order)
ConservedFractions(main, SAGD_df, 4, "SAGD",  all_genes, ct_order)
