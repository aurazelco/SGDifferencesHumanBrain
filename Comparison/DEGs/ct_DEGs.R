# Author: Aura Zelco
# Brief description:
    # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
    # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
    # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
    # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
    # 4. Plots how many genes are found in all age groups, in all but one, etc
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison/DEGs/ct_DEGs_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs/outputs/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/"

# Vectors to save the different sub-groups of DISCO and UCSC
sub_disease <- list.dirs(main_DISCO, full.names = F, recursive = F)
# the first folder "exta_files" is excluded
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco_projs <- c("GSE157827", "GSE174367", "PRJNA544731")
disco <- ImportDataset(main_DISCO, sub_disease, individual_projs = disco_projs)
UCSC <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes")

# disco[[2]] and UCSC[[2]] can be used to manually create unified_annotation, as done below
disco[[2]]
UCSC[[2]]
# manually decided how to combine the sub-celltypes

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

# defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
condition_order <- c("Eze_Nowakowski_integrated_2nd_trimester",
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

# Heatmaps of presence of genes (yes/no) across all ages, for each ct
sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)
PlotCts(main_comparison, sexes, condition_order)

# Heatmaps specifically for the two datasets from the second trimester, all cts
trim_2nd <- CreateConditionDf(c(UCSC[[1]], disco[[1]]), unified_annotation, condition_order[1:2])
PlotAcrossConditions(main_comparison, trim_2nd, "trimester_2nd")


# Count of how genes are shared among ages, for each ct
gene_counts <- CreateCountDfs(sexes)
PlotCountCt(main_comparison, gene_counts)
