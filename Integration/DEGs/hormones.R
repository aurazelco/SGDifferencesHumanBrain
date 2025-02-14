# sources the script containing all functions run here
source("scripts/Integration/DEGs/hormones_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "data/Integration/"

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
hormones <- fromJSON(file="data/Integration/hgv1_hormone_genes.json")
# hormones_source_tgs <- fromJSON(file="Comparison/hgv1_hormone_src_tgt_genes.json")
# hgv1_hormone_src_tgt_genes.json is the same, but the genes are split in source and targets

# Generates a df with all DEGs
sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)

hormones_filt <- hormones[names(which(lapply(hormones, length)>=10))]
df_filt <- CreateHormonesDf(sexes, hormones_filt, groups_order)

hormones_pval <- HormoneEnrichment(df_filt)
write.csv(hormones_pval, paste0(main_comparison, "Hormones/hormone_target_enrichment.csv"))
HmpHormoneEnrichment(main_comparison, hormones_pval, groups_order)
HmpHormoneEnrichment(main_comparison, hormones_pval, groups_order, "Thymosin", "Thymosin")


plot_path <- "data/Integration/Hormones/"
brewer_palette <- brewer.pal(6,"Purples")


hormones_pval_filt <- hormones_pval %>% filter(overlap_hormone_degs > 1)

p <- ggplot(hormones_pval_filt, aes(hormone_id, groups, fill=pval_sign, size=overlap_hormone_degs)) +
  geom_point(color="black", shape=21) +
  facet_grid(ct ~ sex, scales = "free", space = "free") +
  scale_fill_manual(values = c("NS"="white", 
                               "*"=brewer_palette[3],
                               "**"=brewer_palette[4],
                               "***"=brewer_palette[5],
                               "****"=brewer_palette[6]),
                    na.value = "gray") +
  scale_size_binned(range = c(min(hormones_pval_filt$overlap_hormone_degs), 
                                  max(hormones_pval_filt$overlap_hormone_degs))) +
  labs(x="Hormones", y="Groups", size="Number of hormone targets found in SG-DEGs", fill="P-values") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.spacing.x=unit(0.5, "lines"),
        strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
        strip.text.x = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.box = "vertical",
        legend.title = element_text(size=12, face="bold", colour = "black"))



pdf(paste0(plot_path, "Hormones_hypergeom.pdf"), width = 10, height = 14)
print(p)
dev.off()
