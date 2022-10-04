main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/scripts/DEGs/01E_CellMarker_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

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

# NORMAL
PlotCMresults(main, sub_disease[3], 
              "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/", 
              ct_list, data_ct)

# AD
PlotCMresults(main, sub_disease[1], 
              "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/", 
              ct_list, data_ct)
