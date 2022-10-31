main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01B_plot_num_genes_func.R")

####### MAIN

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2

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


# List of diseases sub-folders
sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

# NORMAL
num_df_normal <- IntersectDEG(main, sub_disease[3], pval_thresh, FC_thresh, ct_order)
PlotIntDEGs(main, sub_disease[3], num_df_normal[[1]], num_df_normal[[2]], ct_order)

# AD
num_df_AD <- IntersectDEG(main, sub_disease[1], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[1], num_df_AD[[1]], num_df_AD[[2]], ct_order)

# MS
num_df_MS <- IntersectDEG(main, sub_disease[2], pval_thresh, FC_thresh, ct_order)
PlotIntDEGs(main, sub_disease[2], num_df_MS[[1]], num_df_MS[[2]], ct_order)
