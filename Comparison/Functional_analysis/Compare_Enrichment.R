# Author: Aura Zelco
# Brief description:
  # This script is used to run the enrichment analysis on multiple databases/terms, all using the DEGs from the projects
# Brief procedure:
  # 1. Reads all DEG CSV files from all the different datasets (in this case 2 - DISCO and UCSC)
  # 2. Extract the genes lists for all cts in all conditions
  # 3. For each enrichments, compares F v M in each ct-condition combo and in each sex-ct combo across conditions
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
source("~/Desktop/Lund_MSc/Thesis/scripts/Comparison/Functional_analysis/Compare_Enrichment_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/outputs/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/"

# Need to register for account before running this:
disgenet_api_key <- get_disgenet_api_key(
  email = "aura.zelco@gmail.com", 
  password = "jyqwew-rAbgu3-qyvxuv" )

# Sets the DisGENet API key
Sys.setenv(DISGENET_API_KEY=disgenet_api_key )

# Vectors to save the different sub-groups of DISCO and UCSC
sub_disease <- list.dirs(main_DISCO, full.names = F, recursive = F)
# the first folder "exta_files" is excluded
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco_projs <- c("GSE157827", "GSE174367", "PRJNA544731")
disco <- ImportDataset(main_DISCO, sub_disease, individual_projs = disco_projs, pval = 0.05, FC = 1.2)
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


# Combines all dataframes into one df
sexes <- CreateSexDf(c(disco[[1]], UCSC[[1]][-1]), unified_annotation)

# Saves results, one plot and one CSV for each F v M comparison for each ct-condition combo - GO, KEGG, DO DisGeNET
EnrichFvM(main_comparison, sexes, "GO", "BP")
EnrichFvM(main_comparison, sexes, "KEGG")
EnrichFvM(main_comparison, sexes, "DO")
EnrichFvM(main_comparison, sexes, "DGN")


# Saves results, one plot and one CSV for F and M separately of each ct across conditions - GO, KEGG, DO DisGeNET
EnrichCondition(main_comparison, sexes, "GO", "BP", gene_thresh = 100, condition_ordered = condition_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
EnrichCondition(main_comparison, sexes, "KEGG", gene_thresh = 100, condition_ordered = condition_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
EnrichCondition(main_comparison, sexes, "DO", gene_thresh = 100, condition_ordered = condition_order, rotate_x_axis = T, adj_pval_thresh =  0.01)
EnrichCondition(main_comparison, sexes, "DGN", gene_thresh = 100, condition_ordered = condition_order, rotate_x_axis = T, adj_pval_thresh =  0.01)

# DSigDB - drug db
EnrichOtherDB(main_comparison, sexes, "EnrichR",  "DSigDB", condition_order)
EnrichOtherDBFvM(main_comparison, sexes, "EnrichR",  "DSigDB", condition_order)

# GWAS_Catalog_2019
EnrichOtherDB(main_comparison, sexes, "EnrichR",  "GWAS_Catalog_2019", condition_order)
EnrichOtherDBFvM(main_comparison, sexes, "EnrichR",  "GWAS_Catalog_2019", condition_order)

# DisGeNET (CURATED)
EnrichOtherDB(main_comparison, sexes, "DisGeNET2r",  "DisGeNET (CURATED)", condition_order)
EnrichOtherDBFvM(main_comparison, sexes, "DisGeNET2r",  "DisGeNET (CURATED)", condition_order)

# TANSFAC and JASPAR - TF motifs enrichment
EnrichOtherDB(main_comparison, sexes, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", condition_order)
EnrichOtherDBFvM(main_comparison, sexes, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", condition_order)

# Compare disease-related results
which_comp <- "_comparison_cts"
DO <- ImportDBresults(main_comparison, "DO", which_comp)
DGN <- ImportDBresults(main_comparison, "DGN", which_comp)
DGN_CURATED <- ImportDBresults(main_comparison, "DisGeNET2r_DisGeNET_CURATED", "")
DSigDB <- ImportDBresults(main_comparison, "EnrichR_DSigDB", "")
GWAS <- ImportDBresults(main_comparison, "EnrichR_GWAS_Catalog_2019", "")


# Counts how mauch frequent each term is, adn saves the CSV in the directory
all_dbs <- rbind(DO, DGN, DGN_CURATED, GWAS, DSigDB)
all_dbs <- all_dbs[which(all_dbs$condition!="Eze_Nowakowski_integrated_2nd_trimester"),]
CountDiseases(main_comparison, all_dbs) 

# DGN excluded because too many terms
all_dbs <- rbind(DO, DGN_CURATED, GWAS, DSigDB)
all_dbs <- all_dbs[which(all_dbs$condition!="Eze_Nowakowski_integrated_2nd_trimester"),]
PlotFacetedDB(main_comparison, all_dbs, condition_order )
  
# Cell Enrichment - only adults compared to McKenzie 2018
ct_ref <- as.data.frame(read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/McKenzie_2018_suppl.xlsx",
                    sheet = 'top_human_enrich',
                    skip = 2))

ref_ct_names <- c(
  "ast" = "Astrocytes", 
  "end" = "Endothelial Cells",
  "mic" = "Microglia",
  "neu" = "Neurons",
  "oli" = "Oligodendrocytes"
)

PlotRefCt(main_comparison, sexes, ct_ref, condition_order[8:14], "McKenzie_2018", ref_ct_names)

# Disease enrichment - comparison with Chlamydas 2022
chlamydas <- as.data.frame(read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/Chlamydas_2022.xlsx", skip = 1))
colnames(chlamydas) <- str_replace_all(colnames(chlamydas), " ", "_")
chlamydas <- chlamydas[, c(1,2,4)]
chlamydas <- drop_na(chlamydas)

chl_deg <- CreateDisDf(main_comparison, chlamydas, sexes, "Chlamydas_2022")

PlotDisDeg(main_comparison, chl_deg, "Chlamydas_2022")
