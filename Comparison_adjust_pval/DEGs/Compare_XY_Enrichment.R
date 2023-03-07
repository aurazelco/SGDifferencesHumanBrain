# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the enrichment in the X and Y chromosomes from the DEGs analysis
# Brief procedure:
  # 1. Reads all enrichment csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. 

# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/DEGs/Compare_XY_Enrichment_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison_adjust_pval/"

# Vectors to save the different sub-groups of DISCO and UCSC
# the first folder "exta_files" is excluded
sub_projs <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Common cell type annotation
unified_annotation <- c("CXCL14 IN" = "Interneurons",
                        "EC" = "Endothelial cells",
                        "fibrous astrocyte"  = "Astrocytes",
                        "L2_3 EN" = "Excitatory neurons", 
                        "L4 EN" = "Excitatory neurons",
                        "L5 EN" = "Excitatory neurons",
                        "L5_6 EN" = "Excitatory neurons",
                        "L5b EN" = "Excitatory neurons",
                        "L6 EN" = "Excitatory neurons",                
                        "microglia" = "Microglia", 
                        "Oligodendrocyte" =  "Oligodendrocytes",      
                        "OPC" = "OPCs",                  
                        "PLCH1 L4_5 EN" = "Excitatory neurons", 
                        "protoplasmic astrocyte" = "Astrocytes",
                        "PVALB IN"  = "Interneurons",            
                        "pyramidal neuron"  = "Excitatory neurons",
                        "SST IN" = "Interneurons",   
                        "SV2C IN"  = "Interneurons",   
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
                        "Mesenchymal" = "Mesenchymal",      
                        "Neuroepithelial" =     "Neuroepithelial",
                        "Neuronal" = "Neurons",            
                        "Other"    = "Other",                
                        "Radial Glial"     = "Radial Glia",       
                        "Astrocytes" = "Astrocytes",        
                        "Excitatory neurons"  = "Excitatory neurons",
                        "Interneurons"   = "Interneurons",     
                        "Microglia"  = "Microglia",         
                        "Oligodendrocytes" = "Oligodendrocytes",
                        "OPCs" = "OPCs",            
                        "Unknown" = "Unknown",           
                        "Vascular cells" = "Vascular cells",     
                        "Dorsal progenitors"  = "Dorsal progenitors" ,   
                        "Ventral progenitors" = "Ventral progenitors")
names(unified_annotation) <- tolower(names(unified_annotation))

# Defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
groups_order <- c(
                     "Velmeshev_2022_2nd_trimester",           
                     "Velmeshev_2022_3rd_trimester", 
                     "Velmeshev_2022_0_1_years",                
                     "Velmeshev_2022_1_2_years",            
                     "Velmeshev_2022_2_4_years",  
                     "Velmeshev_2022_10_20_years",      
                     "Velmeshev_2022_Adult",
                     "Normal_GSE157827",              
                     "Normal_GSE174367",               
                     "Normal_PRJNA544731", 
                     "Alzheimer's disease_GSE157827",
                     "Alzheimer's disease_GSE174367",
                     "Multiple Sclerosis_PRJNA544731" 
)

# Imports XY enrichment from all sub-folders
disco <- ImportDataset(main_DISCO, sub_projs, individual_projs = T)
names(disco[[1]]) <- str_replace_all(names(disco[[1]]), "Normal", "Healthy")
UCSC <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes")

# Generates heatmap plot of pvalues and saves it to output folder
PlotEnrichedPvalues(main_comparison, c(UCSC, disco), unified_annotation, groups_order)