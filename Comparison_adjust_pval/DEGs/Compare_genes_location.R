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
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/DEGs/Compare_genes_location_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison_adjust_pval/"

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

# Generates a df with all DEGs
presence <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)

location_df <- readxl::read_xlsx(paste0(main_comparison,"Thul_2017_sm_table_s6.xlsx" ))

genes_loc <- ExtractLocation(presence, location_df, groups_order)

out_path <- paste0(main_comparison, "Genes_location/")
dir.create(out_path, showWarnings = F, recursive = T)

custom_pal <- createPalette(30, c("#010101", "#ff0000"), M=1000)
names(custom_pal) <- c("Uncertain", unique(genes_loc$locations)[1:29])

pdf(paste0(out_path, "Locations.pdf"), height = 15, width = 11)
print(
  ggplot(genes_loc, aes(locations, groups, size=loc_count, color=locations)) +
    geom_point() +
    facet_grid(ct ~ sex, scales = "free", space = "free") +
    scale_colour_manual(values=custom_pal) +
    labs(x="", y="Groups", size="Genes found in location", color="Cellular locations") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.box = "vertical",
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()

all_genes <- do.call(rbind, presence)
all_genes$ct <- gsub("\\..*", "", rownames(all_genes))
all_genes$presence <- str_replace_all(all_genes$presence, c("yes"="Yes", "no"="No"))

mit_genes_ids <- unique(all_genes$gene_id[which(grepl("^MT-", all_genes$gene_id))])
timtom_ids <- c(unique(all_genes$gene_id[which(grepl("^TIMM", all_genes$gene_id))]),
                unique(all_genes$gene_id[which(grepl("^TOMM", all_genes$gene_id))]))
mit_genes_ids <- c("XIST", timtom_ids, mit_genes_ids)

MT_loc <- ExtractLocation(presence, location_df, groups_order, mit_genes_ids)

pdf(paste0(out_path, "MT_genes_locations.pdf"), height = 15, width = 11)
print(
  ggplot(MT_loc, aes(locations, groups, size=loc_count, color=locations)) +
    geom_point() +
    facet_grid(ct ~ sex, scales = "free", space = "free") +
    scale_colour_manual(values=custom_pal) +
    labs(x="", y="Groups", size="Genes found in location", color="Cellular locations") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.box = "vertical",
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()