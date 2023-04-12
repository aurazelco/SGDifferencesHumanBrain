# Author: Aura Zelco
# Brief description:
  # This script is used to run the enrichment analysis on multiple databases/terms, all using the DEGs from the projects
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all groups
  # 3. For each enrichments, compares F v M in each ct-groups combo and in each sex-ct combo across groups
    # clusterProfile:
      # a. GO - Gene Ontology 
      # b. KEGG - Kyoto Encyclopedia of Genes and Genomes
      # c. DO - Disease Ontology
      # d. DGN - DisGeNET
    # enrichR:
      # a. DSigDB
      # b. GWAS_Catalog_2019
      # c. TANSFAC and JASPAR - TF motifs enrichment
    # disgenet2r:
      # a. Curated DisGeNET
    # comparison with McKenzie
    # comparison with Chlamydas
  # 4. Saves the plots and CSV results

# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("~/Desktop/Lund_MSc/Thesis/scripts/Integration/Functional_analysis/Enrichment_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_int_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/"

# Need to register for account before running this:
disgenet_api_key <- get_disgenet_api_key(
  email = "aura.zelco@gmail.com", 
  password = "jyqwew-rAbgu3-qyvxuv" )

# Sets the DisGENet API key
Sys.setenv(DISGENET_API_KEY=disgenet_api_key )

# Vectors to save the different sub-groups of DISCO and UCSC
# the first folder "exta_files" is excluded
sub_projs <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/groups - slightly different file tree structure requires a different approach for UCSC
disco <- ImportDatasets(main_DISCO, sub_projs, UCSC_flag = "no", individual_projs = T)
names(disco[[1]]) <- str_replace_all(names(disco[[1]]), "Normal", "Healthy")
UCSC <- ImportDatasets(main_UCSC, sub_UCSC, UCSC_flag = "yes", individual_projs = F)

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

cts_order <- c(
                      "Excitatory neurons",  
                      "Interneurons",        
                      "Astrocytes",          
                      "Microglia" ,          
                      "Oligodendrocytes",    
                      "OPCs",   
                      "Endothelial cells", 
                      "Vascular cells",      
                      "Dorsal progenitors",  
                      "Ventral progenitors",
                      "Unknown"   
)

# Combines all dataframes into one df
sexes <- CreateSexDf(c(disco[[1]], UCSC[[1]]), unified_annotation)

# Saves results, one plot and one CSV for each F v M comparison for each ct-groups combo - GO, KEGG, DO DisGeNET
#EnrichFvM(main_int_path, sexes, "GO", "BP")
#EnrichFvM(main_int_path, sexes, "KEGG")
#EnrichFvM(main_int_path, sexes, "DO")
#EnrichFvM(main_int_path, sexes, "DGN")


# Saves results, one plot and one CSV for F and M separately of each ct across groups - GO, KEGG, DO DisGeNET
EnrichCts(main_int_path, sexes, "GO", "BP", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
#EnrichCts(main_int_path, sexes, "KEGG", gene_thresh = "no", groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
EnrichCts(main_int_path, sexes, "DO", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
EnrichCts(main_int_path, sexes, "DGN", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01)

# Saves results, one plot and one CSV for F and M separately for each group across cts - GO, KEGG
EnrichGroups(main_int_path, sexes, "GO", "BP", gene_thresh = 100, cts_ordered = cts_order, rotate_x_axis = T, adj_pval_thresh =  0.01)

# DSigDB - drug db
EnrichOtherDB(main_int_path, sexes, "EnrichR",  "DSigDB", groups_order)
#EnrichOtherDBFvM(main_int_path, sexes, "EnrichR",  "DSigDB", groups_order)

# GWAS_Catalog_2019
EnrichOtherDB(main_int_path, sexes, "EnrichR",  "GWAS_Catalog_2019", groups_order)
#EnrichOtherDBFvM(main_int_path, sexes, "EnrichR",  "GWAS_Catalog_2019", groups_order)

# DisGeNET (CURATED)
EnrichOtherDB(main_int_path, sexes, "DisGeNET2r",  "DisGeNET (CURATED)", groups_order)
#EnrichOtherDBFvM(main_int_path, sexes, "DisGeNET2r",  "DisGeNET (CURATED)", groups_order)

# TANSFAC and JASPAR - TF motifs enrichment
EnrichOtherDB(main_int_path, sexes, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", groups_order)
#EnrichOtherDBFvM(main_int_path, sexes, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", groups_order)
EnrichOtherDBGroup(main_int_path, sexes, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", cts_order)


# KEGG - other way instead compareCluster
EnrichOtherDB(main_int_path, sexes, "EnrichR",  "KEGG_2021_Human", groups_order)
#EnrichOtherDBFvM(main_int_path, sexes, "EnrichR",  "KEGG_2021_Human", groups_order)
EnrichOtherDBGroup(main_int_path, sexes, "EnrichR",  "KEGG_2021_Human", cts_order)

# Compare disease-related results
which_comp <- "_comparison_cts"
DO <- ImportDBresults(main_int_path, "DO", which_comp)
DGN <- ImportDBresults(main_int_path, "DGN", which_comp)
DGN_CURATED <- ImportDBresults(main_int_path, "DisGeNET2r_DisGeNET_CURATED", "")
GWAS <- ImportDBresults(main_int_path, "EnrichR_GWAS_Catalog_2019", "")

# Counts how mauch frequent each term is, adn saves the CSV in the directory
all_dbs <- rbind(DO, DGN, DGN_CURATED, GWAS)
all_dbs$dbsx <- str_replace_all(all_dbs$dbsx, c("DisGeNET2r_DisGeNET_CURATED" = "DisGeNET", "EnrichR_GWAS_Catalog_2019" = "GWAS"))
dis_counts_sex <- CountDiseases(main_int_path, all_dbs, which_comp = "sex", min_num = 10) 
dis_counts_ct <- CountDiseases(main_int_path, all_dbs, which_comp = "sex_ct", min_num = 5) 
dis_counts_ct2 <- CountDiseases(main_int_path, all_dbs, which_comp = "sex_ct", min_num = 2) 

dis_counts_dbsx <- CountDiseases(main_int_path, all_dbs,  which_comp = "sex_ct_dbsx", min_num = 5) 


PlotFacetedDBSimplified(main_int_path, dis_counts_sex, which_comp = "sex", min_num = 10)
PlotFacetedDBSimplified(main_int_path, dis_counts_ct, which_comp = "sex_ct", cts_order, min_num = 5)
PlotFacetedDBSimplified(main_int_path, dis_counts_ct2, which_comp = "sex_ct", cts_order, min_num = 2)

# DGN excluded because too many terms
all_dbs <- rbind(DO, DGN_CURATED, GWAS)
PlotFacetedDB(main_int_path, all_dbs, groups_order)

# Drug comparison
drugs <- ImportDBresults(main_int_path, "EnrichR_DSigDB", "")
drugs_counts_sex <- CountDrugs(main_int_path, drugs, which_comp = "sex", min_num = 10) 
drugs_counts_ct <- CountDrugs(main_int_path, drugs, which_comp = "sex_ct", min_num = 5) 

PlotFacetedDrugs(main_int_path, drugs_counts_sex, which_comp = "sex", min_num = 10)
PlotFacetedDrugs(main_int_path, drugs_counts_ct, which_comp = "sex_ct", cts_order, min_num = 5)


  

# TRANSFAC_and_JASPAR_PWMs enrichment
tj_results_cts <- ImportTJPWMs(main_int_path)
CalculateSharedTJPWMs(main_int_path, tj_results_cts, 0.5)
#CalculateSharedTJPWMs(main_int_path, tj_results, 0.75)

tj_results_groups <- ImportTJPWMs(main_int_path, which_folder = "_groups")
CalculateSharedTJPWMs(main_int_path, tj_results_groups, 0.5, which_folder = "_groups")
#CalculateSharedTJPWMs(main_int_path, tj_results_groups, 0.75, which_folder = "_groups")



# Create supplementary files
SaveCPResults(main_int_path, "GO_comparison_cts", "GO_cts")
SaveCPResults(main_int_path, "GO_comparison_groups", "GO_groups")
SaveCPResults(main_int_path, "EnrichR_KEGG_2021_Human", "KEGG_cts")
SaveCPResults(main_int_path, "EnrichR_KEGG_2021_Human_groups", "KEGG_groups")
SaveCPResults(main_int_path, "DGN_comparison_cts", "DGN_cts")
SaveCPResults(main_int_path, "DO_comparison_cts", "DO_cts")
SaveCPResults(main_int_path, "EnrichR_GWAS_Catalog_2019", "GWAS_cts")
SaveCPResults(main_int_path, "DisGeNET2r_DisGeNET_CURATED", "DisGeNET2r_cts")
SaveCPResults(main_int_path, "EnrichR_TRANSFAC_and_JASPAR_PWMs", "TF_cts")
SaveCPResults(main_int_path, "EnrichR_TRANSFAC_and_JASPAR_PWMs_groups", "TF_groups")
SaveCPResults(main_int_path, "EnrichR_DSigDB", "DSigDB_cts")
