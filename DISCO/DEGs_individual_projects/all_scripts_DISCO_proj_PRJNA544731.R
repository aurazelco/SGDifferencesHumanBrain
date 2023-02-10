# 01A_generate_DEGs.R

library(Seurat)

project_id <- "PRJNA544731"

disco_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/"
main <- paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/", project_id, "/")

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))
Idents(disco_filt) <- "proj_sex_disease_ct"

min_cells <- 100

final_groups <- read.csv(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/outputs/final_filt_", min_cells, ".csv"))
final_groups <- final_groups[which(final_groups$proj==project_id), ]
sexes <- unique(final_groups$sex)

for (disease_type in unique(final_groups$disease)) {
  path <- paste0(main, disease_type, "/outputs/01A_DEGs/")
  dir.create(path, showWarnings = FALSE, recursive = T)
  final_filt_disease <- subset(final_groups, disease==disease_type)
  for (ct_type in unique(final_filt_disease$ct)) {
        id1 <- paste(project_id, sexes[1], disease_type, ct_type, sep="_")
        id2 <- paste(project_id, sexes[2], disease_type, ct_type, sep="_")
        if ((id1 %in% unique(final_filt_disease$idents)) & (id2 %in% unique(final_filt_disease$idents))) {
          path_ct <- paste0(path, unique(final_filt_disease[which(final_filt_disease$ct==ct_type), "name_subfolders"]))
          dir.create(path_ct, showWarnings = FALSE)
          deg1 <- FindMarkers(disco_filt, 
                              ident.1 = id1, 
                              ident.2 = id2,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
          write.csv(deg1, paste0(path_ct, "/", sexes[1], ".csv"))
          deg2 <- FindMarkers(disco_filt, 
                              ident.1 = id2, 
                              ident.2 = id1,
                              logfc.threshold = 0.25,
                              min.pct = 0.1,
                              only.pos = TRUE)
          write.csv(deg2, paste0(path_ct, "/", sexes[2],".csv"))
    } 
  }
}

#write.csv(final_groups, paste0(main, "final_filt_", min_cells, ".csv"),row.names = F)

rm(list=ls())

####### GENERAL VARIABLES

project_id <- "PRJNA544731"
main <- paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/", project_id, "/")

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

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/01B_plot_num_genes_func.R")

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

# NORMAL
CountDEG(main, sub_disease[2], pval_thresh, FC_thresh, ct_order)

# MS
CountDEG(main, sub_disease[1], pval_thresh, FC_thresh, ct_order)

####### 01C_num_chr.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects//01C_num_chr_func.R")

# NORMAL
chr_normal <- ProcessCt(main, sub_disease[2])
ExtractSharedGenes(main, sub_disease[2], chr_normal)
PlotGeneralHeatmap(main, sub_disease[2], chr_normal, ct_order)
PlotSexHmp(main, sub_disease[2], chr_normal, ct_order)

# MS
chr_ms <- ProcessCt(main, sub_disease[1])
ExtractSharedGenes(main, sub_disease[1], chr_ms)
PlotGeneralHeatmap(main, sub_disease[1], chr_ms, ct_order)
PlotSexHmp(main, sub_disease[1], chr_ms, ct_order)

####### 01D_Xpar1,2.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects//01D_Xpar1,2_func.R")

Xpar1 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

# NORMAL
XparCt(main, sub_disease[2], Xpar1_list, Xpar2_list, ct_order)

# MS
XparCt(main, sub_disease[1], Xpar1_list, Xpar2_list, ct_order)


####### 01E_CellMarker.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/01E_CellMarker_func.R")

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

extra_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/"

# NORMAL
PlotCMresults(main, sub_disease[2], 
              extra_path, 
              ct_list, data_ct, ct_order)

# MS
PlotCMresults(main, sub_disease[1], 
              extra_path, 
              ct_list, data_ct, ct_order)


####### 02A_HyperGeom.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/01C_num_chr_func.R")
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/02A_HyperGeom_func.R")

# as used in 02A_HyperGeom
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

# NORMAL
SexChr2(main, sub_disease[2], tot_genes, X_chr_genes, Y_chr_genes)
PlotNumChr(main, sub_disease[2], num_chr_genes, ct_order, T)

# MS
SexChr2(main, sub_disease[1], tot_genes, X_chr_genes, Y_chr_genes)
PlotNumChr(main, sub_disease[1], num_chr_genes, ct_order, T)

####### 02B_ARE_ERE.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/02B_ARE_ERE_func.R")

ARE <- read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

# NORMAL
AnalysisARE_ERE(main, sub_disease[2], ARE, EREgene, ct_order)

# MS
AnalysisARE_ERE(main, sub_disease[1], ARE, EREgene, ct_order)

####### 02C_Conservation.R

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs_individual_projects/02C_Conservation_func.R")

# CONSERVATION ACROSS PRIMATES
conserved <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/extra_files/mart_export.txt",
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

SAGD <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/extra_files/Sexassociatedgene_Padj0.05_PMID30380119.csv")
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

write.csv(SAGD_df, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/extra_files/SAGD_filt.csv")

SAGD_df <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/extra_files/SAGD_filt.csv")
SAGD_df[,1] <- NULL

# ENSEMBL 

ensembl_mat <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/extra_files/binary_mat_all_species.csv")
names(ensembl_mat)[names(ensembl_mat) == 'X'] <- "gene_name"

# Read all DEGs from healthy patients, regardless of the sex

# all genes commonly expressed in the cts
all_genes <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj/tot_genes_ct.csv")
all_genes$X <- NULL
all_genes <- subset(all_genes, proj==project_id)
all_genes$proj <- NULL
col_factors <- c("disease", "sex", "ct")
all_genes[col_factors] <- lapply(all_genes[col_factors], as.factor) 

# NORMAL
ConservedFractions(main, sub_disease[2], conserved, 4, "Primates", all_genes, ct_order)
ConservedFractions(main, sub_disease[2], SAGD_df, 4, "SAGD",  all_genes, ct_order)
ConservedFractions(main, sub_disease[2], ensembl_mat, 4, "ENSEMBL",  all_genes, ct_order)
ConservedFractions(main, sub_disease[2], ensembl_mat, 100, "ENSEMBL",  all_genes, ct_order)

# MS
ConservedFractions(main, sub_disease[1], conserved, 4, "Primates",  all_genes, ct_order)
ConservedFractions(main, sub_disease[1], SAGD_df, 4, "SAGD",  all_genes, ct_order)
ConservedFractions(main, sub_disease[1], ensembl_mat, 4, "ENSEMBL",  all_genes, ct_order)
ConservedFractions(main, sub_disease[1], ensembl_mat, 100, "ENSEMBL",  all_genes, ct_order)
