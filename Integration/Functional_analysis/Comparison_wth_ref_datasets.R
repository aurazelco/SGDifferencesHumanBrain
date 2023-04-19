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
source("~/Desktop/Lund_MSc/Thesis/scripts/Integration/Functional_analysis/Comparison_wth_ref_datasets_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_int_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/"

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


# Cell Enrichment - only adults compared to McKenzie 2018
ct_ref <- as.data.frame(read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/McKenzie_2018_suppl.xlsx",
                                  sheet = 'top_human_enrich',
                                  skip = 2))

ref_ct_names <- c(
  "ast" = "Astrocytes", 
  "end" = "Endothelial Cells",
  "mic" = "Microglia",
  "neu" = "Neurons",
  "oli" = "Oligodendrocytes"
)

PlotRefCt(main_int_path, sexes, ct_ref, groups_order[7:13], "McKenzie_2018", ref_ct_names, facets_align = "v")
PlotRefCt(main_int_path, sexes, ct_ref, groups_order[7:13], "McKenzie_2018", ref_ct_names, facets_align = "h")


# Disease enrichment - comparison with Chlamydas 2022
chlamydas <- as.data.frame(read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/Chlamydas_2022.xlsx", skip = 1))
colnames(chlamydas) <- str_replace_all(colnames(chlamydas), " ", "_")
chlamydas <- chlamydas[, c(1,2,4)]
chlamydas <- drop_na(chlamydas)

chl_deg <- CreateDisDf(main_int_path, chlamydas, sexes, "Chlamydas_2022")

PlotDisDeg(main_int_path, chl_deg, "Chlamydas_2022", groups_order)

# SFARI autism db
sfari <- read.csv(paste0(main_int_path, "SFARI-Gene_genes.csv"))
count_sfari <- CountSFARI(main_int_path, sexes, sfari, groups_order)
sfari_chr_gene_counts <- c("Autosome" = nrow(sfari[!is.na(as.numeric(sfari$chromosome)),]),
                           "X" = nrow(sfari[which(sfari$chromosome=="X"), ]),      
                           "X,Y" = nrow(sfari[which(sfari$chromosome=="X,Y"), ]),  
                           "Y" = nrow(sfari[which(sfari$chromosome=="Y"), ])
                           
)

PlotSFARI(main_int_path, count_sfari)
PlotSFARI(main_int_path, count_sfari, which_comp = "chr")

tot_genes <- 20000
enriched_sfari_chr <- HyperGeomSFARI(main_int_path, count_sfari, tot_genes, sfari_chr_gene_counts, chr_comp = T)
PlotEnrichedPvalues(main_int_path, enriched_sfari_chr, groups_order, cts_order, chr_comp = T)

enriched_sfari_tot <- HyperGeomSFARI(main_int_path, count_sfari, tot_genes, sfari_chr_gene_counts, chr_comp = F)
PlotEnrichedPvalues(main_int_path, enriched_sfari_tot, groups_order, cts_order, chr_comp = F)

