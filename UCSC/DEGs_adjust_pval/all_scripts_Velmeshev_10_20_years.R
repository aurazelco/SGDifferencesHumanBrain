# 01A_generate_DEGs.R

library(Seurat)
library(stringr)

main <- "UCSC/DEGs_adjust_pval/"

dir.create(main, recursive = T, showWarnings = F)

velm_10_20_years <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_10_20_years.rds")
#velm_10_20_years@meta.data$sex_ct <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$cluster_final, sep="_")
Idents(velm_10_20_years) <- "sex_ct"

min_cells <- 100

main_deg <- paste0(main, "Velmeshev_2022_10_20_years/outputs/")
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
    deg1 <- FindMarkers(velm_10_20_years, 
                              ident.1 = id1, 
                              ident.2 = id2,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
    write.csv(deg1, paste0(path_ct, names(sexes[1]), ".csv"))
    deg2 <- FindMarkers(velm_10_20_years, 
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

main_local <- "UCSC/DEGs_adjust_pval/Velmeshev_2022_10_20_years/outputs/"

ct_order <- c(
  "Dorsal progenitors",
  "Ventral progenitors",
  "Excitatory neurons",    
  "Interneurons",                      
  "OPCs",  
  "Oligodendrocytes",
  "Astrocytes",
  "Microglia",   
  "Vascular cells",  
  "Unknown"
)

####### 01B_plot_num_genes.R

source("scripts/UCSC/DEGs_adjust_pval/01B_plot_num_genes_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2
min_num_thresh <- 10

# NORMAL
CountDEG(main_local, pval_thresh, FC_thresh, ct_order, min_num_thresh)

####### 01C_num_chr.R

source("scripts/UCSC/DEGs_adjust_pval/01C_num_chr_func.R")

chr_10_20_years <- ProcessCt(main_local)
ExtractSharedGenes(main_local, chr_10_20_years)
PlotGeneralHeatmap(main_local, chr_10_20_years, ct_order, "Velmeshev_10_20_years")
PlotSexHmp(main_local, chr_10_20_years, ct_order)

####### 01D_Xpar1,2.R

source("scripts/UCSC/DEGs_adjust_pval/01D_Xpar1,2_func.R")

Xpar1 <- read.csv("UCSC/DEGs_adjust_pval/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("UCSC/DEGs_adjust_pval/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

XparCt(main_local, Xpar1_list, Xpar2_list, ct_order)


####### 02A_HyperGeom.R

source("scripts/UCSC/DEGs_adjust_pval/01C_num_chr_func.R")
source("scripts/UCSC/DEGs_adjust_pval/02A_HyperGeom_func.R")

# as used in 02A_HyperGeom
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

SexChr2(main_local, tot_genes, X_chr_genes, Y_chr_genes)
PlotNumChr(main_local, num_chr_genes, ct_order, T)



####### 02B_ARE_ERE.R

source("scripts/UCSC/DEGs_adjust_pval/02B_ARE_ERE_func.R")

ARE <- read_excel("UCSC/DEGs_adjust_pval/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("UCSC/DEGs_adjust_pval/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

AnalysisARE_ERE(main_local, ARE, EREgene, ct_order)

####### 02C_Conservation.R

source("scripts/UCSC/DEGs_adjust_pval/02C_Conservation_func.R")

# CONSERVATION ACROSS PRIMATES
conserved <- read.csv("UCSC/DEGs_adjust_pval/extra_files/mart_export.txt",
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
SAGD <- read.csv("UCSC/DEGs_adjust_pval/extra_files/Sexassociatedgene_Padj0.05_PMID30380119.csv")
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
names(SAGD_df)[names(SAGD_df) == 'gene'] <- 'gene_name'

write.csv(SAGD_df, "UCSC/DEGs_adjust_pval/extra_files/SAGD_filt.csv")

SAGD_df <- read.csv("UCSC/DEGs_adjust_pval/extra_files/SAGD_filt.csv")
SAGD_df$X <- NULL

# ENSEMBL 

ensembl_mat <- read.csv("UCSC/DEGs_adjust_pval/extra_files/binary_mat_all_species.csv")
names(ensembl_mat)[names(ensembl_mat) == 'X'] <- "gene_name"

# all genes commonly expressed in the cts
all_genes <- read.csv("UCSC/DEGs_adjust_pval/Velmeshev_2022_10_20_years/tot_genes_ct_Velmeshev_2022_10_20_years.csv")
all_genes$X <- NULL
all_genes$sex <- str_replace_all(all_genes$sex, c("Female"="F", "Male"="M"))
col_factors <- c("sex", "ct")
all_genes[col_factors] <- lapply(all_genes[col_factors], as.factor) 


ConservedFractions(main_local, conserved, 4, "Primates", all_genes, ct_order)
ConservedFractions(main_local, SAGD_df, 4, "SAGD",  all_genes, ct_order)
ConservedFractions(main_local, ensembl_mat, 4, "ENSEMBL",  all_genes, ct_order)
ConservedFractions(main_local, ensembl_mat, 100, "ENSEMBL",  all_genes, ct_order)
