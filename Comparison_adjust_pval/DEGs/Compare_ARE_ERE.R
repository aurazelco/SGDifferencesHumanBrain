# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the AREs and EREs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
  # 1. Reads all ARE csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
  # 3. Plots the percentages of ARE sites in each ct across conditions, separated by sex

# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/DEGs/Compare_ARE_ERE_func.R")

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
                     "Healthy_GSE157827",              
                     "Healthy_GSE174367",               
                     "Healthy_PRJNA544731", 
                     "Alzheimer's disease_GSE157827",
                     "Alzheimer's disease_GSE174367",
                     "Multiple Sclerosis_PRJNA544731" 
)

# Imports ARE from all sub-folders
disco_ARE <- ImportDataset(main_DISCO, sub_projs, individual_projs = T, ARE_ERE="ARE")
names(disco_ARE) <- str_replace_all(names(disco_ARE), "Normal", "Healthy")
UCSC_ARE <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes", ARE_ERE="ARE")

# Combines them in one dataframe (summing the common annotation) and plots the results
ARE <- CreateAREDf(c(disco_ARE, UCSC_ARE), unified_annotation)
PlotARE(main_comparison, ARE, groups_order)
# Facets all ARE plots
PlotFacetedARE(main_comparison, ARE, groups_order)

# Plot simplified sites facet
ARE_simpl <- CreateAREDf(c(disco_ARE, UCSC_ARE), unified_annotation, "yes")
PlotFacetedARE(main_comparison, ARE_simpl, groups_order, simpl = "yes")



# Imports ERE from all sub-folders
disco_ERE <- ImportDataset(main_DISCO, sub_projs, individual_projs = T, ARE_ERE="ERE")
names(disco_ERE) <- str_replace_all(names(disco_ERE), "Normal", "Healthy")

UCSC_ERE <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes", ARE_ERE="ERE")

# Combines them in one dataframe (summing the common annotation) and plots the results
ERE <- CreateEREDf(c(disco_ERE, UCSC_ERE), unified_annotation)
# Facets all ERE plots
PlotFacetedERE(main_comparison, ERE, groups_order)
