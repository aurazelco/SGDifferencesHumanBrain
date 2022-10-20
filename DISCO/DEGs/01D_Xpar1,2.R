main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/01D_Xpar1,2_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

Xpar1 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar1.csv",
                  skip = 1)
Xpar1_list <- Xpar1$Approved.symbol
Xpar2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Xpar2.csv",
                  skip = 1)
Xpar2_list <- Xpar2$Approved.symbol

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
XparCt(main, sub_disease[3], Xpar1_list, Xpar2_list, ct_order)


#AD
XparCt(main, sub_disease[1], Xpar1_list, Xpar2_list, ct_order)

