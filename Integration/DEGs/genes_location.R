# Author: Aura Zelco
# Brief description:
    # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
    # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
    # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
    # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
    # 4. Plots how many genes are found in all age groups, in all but one, etc
    # 5. Plots the total number of DEGs per ct across all conditions 
    # 6. Plots the number of overlapping genes between a specific condition and all others, divided by ct and sex
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Integration/DEGs/genes_location_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Integration/"

# Vectors to save the different sub-groups of DISCO and UCSC
sub_disco <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
# the first folder "exta_files" is excluded
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco <- ImportDatasets(main_DISCO, sub_disco, UCSC_flag = "no", individual_projs = T)
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

# Generates a df with all DEGs
presence <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)

location_df <- readxl::read_xlsx(paste0(main_comparison,"Thul_2017_sm_table_s6.xlsx" ))

genes_loc <- ExtractLocation(presence, location_df, groups_order)

PlotLocation(main_comparison, genes_loc, plot_title = "Cellular_compartments")

loc_counts <- c(colSums(location_df[which(location_df$Reliability!="Uncertain"), c(4:32)]),
                "Uncertain"=sum(colSums(location_df[which(location_df$Reliability=="Uncertain"), c(4:32)])))

tot_genes <- 20000

enriched_loc <- HyperGeomLocation(main_comparison, genes_loc, tot_genes, loc_counts)


location_results <- merge(genes_loc, enriched_loc, by=c("groups", "ct", "sex", "locations"))


PlotEnrichedPvalues(main_comparison, location_results, groups_order, cts_order, plot_type="dot")
PlotEnrichedPvalues(main_comparison, location_results, groups_order, cts_order, plot_type="hmp")

# Plot location of MT genes
all_genes <- do.call(rbind, presence)
all_genes$ct <- gsub("\\..*", "", rownames(all_genes))
all_genes$presence <- str_replace_all(all_genes$presence, c("yes"="Yes", "no"="No"))

mit_genes_ids <- unique(all_genes$gene_id[which(grepl("^MT-", all_genes$gene_id))])
timtom_ids <- c(unique(all_genes$gene_id[which(grepl("^TIMM", all_genes$gene_id))]),
                unique(all_genes$gene_id[which(grepl("^TOMM", all_genes$gene_id))]))
mit_genes_ids <- c("XIST", timtom_ids, mit_genes_ids)

MT_loc <- ExtractLocation(presence, location_df, groups_order, mit_genes_ids)
PlotLocation(main_comparison, MT_loc, "MT_genes_locations")
