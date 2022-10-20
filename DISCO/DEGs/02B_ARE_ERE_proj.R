main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02B_ARE_ERE_proj_func.R")

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

ARE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

# List of diseases sub-folders
sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

# NORMAL
AnalysisARE_ERE(main, sub_disease[3], pval_thresh, FC_thresh, ARE, EREgene, ct_order)

# AD
AnalysisARE_ERE(main, sub_disease[1], pval_thresh, FC_thresh, ARE, EREgene, ct_order)









