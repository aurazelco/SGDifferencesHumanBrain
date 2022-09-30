main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01C_num_chr_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

# as used in 02A_Fisher
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000

num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

# NORMAL
normal <- sub_disease[3]
chr_normal <- ProcessCt(main, normal)
PlotSexHmp(main, normal, chr_normal)
PlotNumChr(main, normal, num_chr_genes, T)

#AD
ad <- sub_disease[1]
chr_ad <- ProcessCt(main, ad)
PlotSexHmp(main, ad, chr_ad)
PlotNumChr(main, ad, num_chr_genes, T)


