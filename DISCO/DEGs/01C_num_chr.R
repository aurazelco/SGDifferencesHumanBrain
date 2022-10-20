main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01C_num_chr_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

# as used in 02A_Fisher
X_chr_genes <- 1848
Y_chr_genes <- 431
tot_genes <- 20000
num_chr_genes <- list("X" = X_chr_genes, "Y" = Y_chr_genes, "Autosome" = (tot_genes - X_chr_genes - Y_chr_genes))

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

# NORMAL
chr_normal <- ProcessCt(main, sub_disease[3])
PlotGeneralHeatmap(main, sub_disease[3], chr_normal, ct_order)
PlotSexHmp(main, sub_disease[3], chr_normal, ct_order)
PlotNumChr(main, sub_disease[3], num_chr_genes, T, ct_order)



#AD
chr_ad <- ProcessCt(main, sub_disease[1])
PlotGeneralHeatmap(main, sub_disease[1], chr_ad, ct_order)
PlotSexHmp(main, sub_disease[1], chr_ad, ct_order)
PlotNumChr(main, sub_disease[1], num_chr_genes, T, ct_order)


