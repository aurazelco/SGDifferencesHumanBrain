# Author: Aura Zelco
# Brief description:
  # This script is used to plot the DEGs from UCSC 2nd trimester and O'Brien as Venn diagrams
# Brief procedure:
  # 1. Reads all DEG csv files from all the different datasets
  # 2. Imports the reference dataset - in this case O'Brien et al. 2019
  # 3. Plots the resulting Venn diagrams, one for each sex
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Integration/DEGs/DEGs_2nd_trimester_func.R")

# sets the directories where to find the DEG csv files
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_2nd_trimester/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_int_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/"

# the 2nd_trimester folders are selected
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)
sub_UCSC <- sub_UCSC[which(grepl("2nd_trimester", sub_UCSC))]

pval_thresh <- 0.05
FC_thresh <- 1.2

# Import all the CSVs from the different sub-folders
UCSC <- ImportDataset(main_UCSC, sub_UCSC, pval = pval_thresh, FC = FC_thresh)

# Import the reference from the bioRXiv paper O'Brien et al. 2019
obrien <- read_xlsx(paste0(main_int_path, "O'Brien_bioRXiv_2019_suppl.xlsx"), sheet = 2, skip = 1)
colnames(obrien) <- str_replace_all(colnames(obrien), " ", "_")

# Extract the genes from each list and plots the Venn diagrams, one per sex, plus the percentages of overlap
#Venn2ndTrimSex(main_int_path, obrien, n_col = 2, UCSC, "O'Brien",to_remove = "_2nd_trimester")

# Using Eze as a reference
comps <- c( "O'Brien", "Eze_Nowakowski_integrated_F", "Eze_Nowakowski_integrated_M", "Velmeshev_2022_F", "Velmeshev_2022_M" )
name_comps <- c("O'Brien", "Eze-Nowakowski F", "Eze-Nowakowski M", "2nd trimester F", "2nd trimester M")

Venn2ndTrim(main_int_path, obrien, list_degs = UCSC, ref_name =  "O'Brien", groups_to_compare = comps[c(1, 4, 5)], plot_title = "Velmeshev_vs_O'Brien", to_remove = "_2nd_trimester", cat_names = name_comps[c(1, 4, 5)])
Venn2ndTrim(main_int_path, obrien, list_degs = UCSC, ref_name =  "O'Brien", groups_to_compare = comps[c(2, 4)], plot_title = "Velmeshev_vs_Eze-Nowakowski_F", to_remove = "_2nd_trimester", cat_names = name_comps[c(2, 4)])
Venn2ndTrim(main_int_path, obrien, list_degs = UCSC, ref_name =  "O'Brien", groups_to_compare = comps[c(3, 5)], plot_title = "Velmeshev_vs_Eze-Nowakowski_M", to_remove = "_2nd_trimester", cat_names = name_comps[c(3,5)])
Venn2ndTrim(main_int_path, obrien, list_degs = UCSC, ref_name =  "O'Brien", groups_to_compare = comps, plot_title = "All", to_remove = "_2nd_trimester", cat_names = name_comps)

