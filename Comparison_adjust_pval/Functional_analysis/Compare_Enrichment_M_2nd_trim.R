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
source("~/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/Functional_analysis/Compare_Enrichment_M_2nd_trim_func.R")

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison_adjust_pval/"

# Need to register for account before running this:
disgenet_api_key <- get_disgenet_api_key(
  email = "aura.zelco@gmail.com", 
  password = "jyqwew-rAbgu3-qyvxuv" )

# Sets the DisGENet API key
Sys.setenv(DISGENET_API_KEY=disgenet_api_key )

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
  "Microglia",          
  "OPCs",                
  "Vascular cells",  
  "Dorsal progenitors",  
  "Ventral progenitors"
)

# Import shared DEGs for M second trimester
M_shared <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/Velmeshev_2022_2nd_trimester/outputs/01C_num_chr/M_shared_genes.csv")
M_shared$X <- NULL
M_shared_ls <- split(M_shared[, c(1,2)], M_shared$ct, drop = T) 
M_shared_ls <- lapply(M_shared_ls, function(x) { x["ct"] <- NULL; x })
M_shared_ls <- lapply(M_shared_ls, function(x) { x <- x[["gene"]] })
M_shared_ls <- M_shared_ls[cts_order[which(cts_order %in% names(M_shared_ls))]]

Enrich2ndTrim(main_comparison, M_shared_ls, "GO", "BP", gene_thresh = 100, rotate_x_axis = T, adj_pval_thresh =  0.01)
Enrich2ndTrim(main_comparison, M_shared_ls, "DO", gene_thresh = 100, rotate_x_axis = T, adj_pval_thresh =  0.01)
Enrich2ndTrim(main_comparison, M_shared_ls, "DGN", gene_thresh = 100, rotate_x_axis = T, adj_pval_thresh =  0.01)

# DSigDB - drug db
EnrichOtherDBGroup2ndTrim(main_comparison, M_shared_ls, "EnrichR",  "DSigDB", cts_order)

# GWAS_Catalog_2019
EnrichOtherDBGroup2ndTrim(main_comparison, M_shared_ls, "EnrichR",  "GWAS_Catalog_2019", cts_order)

# DisGeNET (CURATED)
EnrichOtherDBGroup2ndTrim(main_comparison, M_shared_ls, "DisGeNET2r",  "DisGeNET (CURATED)", cts_order)

# TANSFAC and JASPAR - TF motifs enrichment
EnrichOtherDBGroup2ndTrim(main_comparison, M_shared_ls, "EnrichR",  "TRANSFAC_and_JASPAR_PWMs", cts_order)

# KEGG - other way instead compareCluster
EnrichOtherDBGroup2ndTrim(main_comparison, M_shared_ls, "EnrichR",  "KEGG_2021_Human", cts_order)
