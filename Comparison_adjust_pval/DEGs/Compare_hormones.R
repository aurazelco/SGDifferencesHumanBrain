# sources the script containing all functions run here
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/DEGs/Compare_hormones_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison_adjust_pval/"

# Vectors to save the different sub-groups of DISCO and UCSC
# the first folder "exta_files" is excluded
sub_projs <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
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

# Imports the hormone file
hormones <- fromJSON(file="/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison_adjust_pval/hgv1_hormone_genes.json")
# hormones_source_tgs <- fromJSON(file="/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/hgv1_hormone_src_tgt_genes.json")
# hgv1_hormone_src_tgt_genes.json is the same, but the genes are split in source and targets

# Generates a df with all DEGs
sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)

hormones_filt <- hormones[names(which(lapply(hormones, length)>=10))]
df_filt <- CreateHormonesDf(sexes, hormones_filt, groups_order)

PlotHormonesRes(main_comparison, df_filt, groups_order, "abs")
PlotHormonesRes(main_comparison, df_filt, groups_order, "perc_degs")
PlotHormonesRes(main_comparison, df_filt, groups_order, "perc_hormones")

PlotHormonesResFaceted(main_comparison, df_filt, groups_order, "abs")
PlotHormonesResFaceted(main_comparison, df_filt, groups_order, "perc_degs")
PlotHormonesResFaceted(main_comparison, df_filt, groups_order, "perc_hormones")

hormones_pval <- HormoneEnrichment(df_filt)
HmpHormoneEnrichment(main_comparison, hormones_pval, groups_order)
HmpHormoneEnrichment(main_comparison, hormones_pval, groups_order, "Thymosin", "Thymosin")


intersect(
  tolower(unique(sexes[which(sexes$sex=="F" & sexes$common_annot=="Microglia" & sexes$condition=="Alzheimer's disease_GSE157827"), "gene_id"])),
  hormones_filt$testosterone
)
# "dusp1" "calr" 
intersect(
  tolower(unique(sexes[which(sexes$sex=="F" & sexes$common_annot=="Microglia" & sexes$condition=="Alzheimer's disease_GSE174367"), "gene_id"])),
  hormones_filt$testosterone
)
#  "spp1"