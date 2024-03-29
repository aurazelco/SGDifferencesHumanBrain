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
source("scripts/Integration/Functional_analysis/Comparison_wth_ref_datasets_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_int_path <- "Integration/"

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
ct_ref <- as.data.frame(read_xlsx(paste0(main_int_path, "McKenzie_2018_suppl.xlsx"),
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
chlamydas <- as.data.frame(read_xlsx(paste0(main_int_path, "Chlamydas_2022.xlsx"), skip = 1))
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

# Comparison withOliva et al 2020
tab_names <- excel_sheets(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"))
brain_tabs <- c(tab_names[which(grepl("^BRN", tab_names))], "PTTARY")

oliva <- lapply(brain_tabs, function(x) read_excel(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"), sheet = x))
names(oliva) <- brain_tabs
oliva <- do.call(rbind, oliva)
oliva$region <- gsub("\\..*", "", rownames(oliva))
rownames(oliva) <- NULL
oliva <- as.data.frame(oliva)
oliva$chr_simplified <- str_replace_all(oliva$chr, c("chrX"="X", "chrY"="Y", "chr\\d+"="Autosome"))
oliva <- oliva[which(!is.na(oliva$HUGO_gene_id)), ]

oliva_num <- c("X"=length(unique(oliva[which(oliva$chr_simplified=="X"), "HUGO_gene_id"])),
               "Autosome" = length(unique(oliva[which(oliva$chr_simplified=="Autosome"), "HUGO_gene_id"])))
tot_genes <- 20000

oliva_num_reg <- c("BRNAMY" = length(oliva[which(oliva$region=="BRNAMY"), "HUGO_gene_id"]), 
                  "BRNACC"  = length(oliva[which(oliva$region=="BRNACC"), "HUGO_gene_id"]),
                  "BRNCDT"  = length(oliva[which(oliva$region=="BRNCDT"), "HUGO_gene_id"]),
                  "BRNCHB"  = length(oliva[which(oliva$region=="BRNCHB"), "HUGO_gene_id"]),
                  "BRNCHA"  = length(oliva[which(oliva$region=="BRNCHA"), "HUGO_gene_id"]),
                  "BRNCTXA" = length(oliva[which(oliva$region=="BRNCTXA"), "HUGO_gene_id"]),
                  "BRNCTXB" = length(oliva[which(oliva$region=="BRNCTXB"), "HUGO_gene_id"]),
                  "BRNHPP"  = length(oliva[which(oliva$region=="BRNHPP"), "HUGO_gene_id"]),
                  "BRNHPT"  = length(oliva[which(oliva$region=="BRNHPT"), "HUGO_gene_id"]),
                  "BRNNCC"  = length(oliva[which(oliva$region== "BRNNCC"), "HUGO_gene_id"]),
                  "BRNPTM"  = length(oliva[which(oliva$region=="BRNPTM"), "HUGO_gene_id"]),
                  "BRNSPC"  = length(oliva[which(oliva$region=="BRNSPC"), "HUGO_gene_id"]),
                  "BRNSNG" = length(oliva[which(oliva$region=="BRNSNG"), "HUGO_gene_id"]),
                  "PTTARY" = length(oliva[which(oliva$region=="PTTARY"), "HUGO_gene_id"])
                  )


# On all regions together
count_oliva_all <- CountOliva(main_int_path, sexes, oliva, groups_order)
count_oliva_all$oliva_perc <- count_oliva_all$oliva_count * 100 / count_oliva_all$tot_degs_count

oliva_pval <- HyperGeomOliva(main_int_path, count_oliva_all, tot_genes, oliva_num, chr_comp = T)
oliva_pval_all <- HyperGeomOliva(main_int_path, count_oliva_all, tot_genes, oliva_num, chr_comp = F)

PlotEnrichedPvaluesOliva(main_int_path, oliva_pval, groups_order, cts_order, chr_comp = T)
PlotEnrichedPvaluesOliva(main_int_path, oliva_pval_all, groups_order, cts_order, chr_comp = F)


# separated by region
count_oliva_reg <- CountOliva(main_int_path, sexes, oliva, groups_order, reg_split = T)
count_oliva_reg$oliva_perc <- count_oliva_reg$oliva_count * 100 / count_oliva_reg$tot_degs_count

oliva_pval_reg <- HyperGeomOlivaReg(main_int_path, count_oliva_reg, tot_genes, oliva_num_reg)

PlotEnrichedPvaluesOlivaReg(main_int_path, oliva_pval_reg, groups_order, cts_order)
