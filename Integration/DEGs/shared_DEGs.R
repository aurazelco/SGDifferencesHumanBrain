# Author: Aura Zelco
# Brief description:
  # This script is used for comparing the shared DEGs across the cts from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
  # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Manually combines the annotations to be able to compare at a general level the different groups
  # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each chromosome and sex
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("scripts/Integration/DEGs/shared_DEGs_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "Integration/"

# Vectors to save the different sub-groups of DISCO and UCSC
# the first folder "extra_files" is excluded
sub_projs <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco <- ImportDataset(main_DISCO, sub_projs, individual_projs = T)
names(disco) <- str_replace_all(names(disco), "Normal", "Healthy")
UCSC <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes")

# defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
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

# Re-organizes the input files into one
all_shared <- CreateDf(c(UCSC, disco))

# Plots the DEGs according to their chromosome, and if they are found or not in each group
PlotSharedDEGs(main_comparison, all_shared, "Autosome", groups_order)
PlotSharedDEGs(main_comparison, all_shared, "X", groups_order)
PlotSharedDEGs(main_comparison, all_shared, "Y", groups_order)