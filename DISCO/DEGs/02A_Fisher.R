main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02A_Fisher_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

# NORMAL
normal <- sub_disease[3]
SexChr(main, normal, tot_genes, X_chr_genes, Y_chr_genes)
SexChr2(main, normal, tot_genes, X_chr_genes, Y_chr_genes)




#AD
ad <- sub_disease[1]
#SexChr(main, ad, X_chr_genes, Y_chr_genes)
SexChr2(main, ad, tot_genes, X_chr_genes, Y_chr_genes)

