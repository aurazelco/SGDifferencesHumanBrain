main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/scripts/DEGs/02C_Conservation_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

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



# NORMAL
normal_all_deg <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/normal_all_DEGs.csv")
ConservedFractions(main, sub_disease[3], conserved, 4, "Primates", normal_all_deg)
ConservedFractions(main, sub_disease[3], SAGD_df, 4, "SAGD", normal_all_deg)


#AD
ad_all_deg <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/ad_all_DEGs.csv")
ConservedFractions(main, sub_disease[1], conserved, 4, "Primates", ad_all_deg)
ConservedFractions(main, sub_disease[1], SAGD_df, 4, "SAGD", ad_all_deg)


