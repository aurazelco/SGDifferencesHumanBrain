main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/scripts/DEGs/01B_plot_num_genes_func.R")

####### MAIN

# QC parameters
pval_thresh <- 0.05
FC_thresh <- 1.2


# List of diseases sub-folders
sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

# NORMAL
num_df_normal <- IntersectDEG(main, sub_disease[3], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[3], num_df_normal[[1]], num_df_normal[[2]])



# AD
num_df_AD <- IntersectDEG(main, sub_disease[1], pval_thresh, FC_thresh)
PlotIntDEGs(main, sub_disease[1], num_df_AD[[1]], num_df_AD[[2]])




