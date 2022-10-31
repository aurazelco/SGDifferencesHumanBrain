main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02A_Fisher_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

# NORMAL
#SexChr(main, sub_disease[3], tot_genes, X_chr_genes, Y_chr_genes)
SexChr2(main, sub_disease[3], tot_genes, X_chr_genes, Y_chr_genes)

# AD
SexChr2(main, sub_disease[1], tot_genes, X_chr_genes, Y_chr_genes)

# MS
SexChr2(main, sub_disease[2], tot_genes, X_chr_genes, Y_chr_genes)
