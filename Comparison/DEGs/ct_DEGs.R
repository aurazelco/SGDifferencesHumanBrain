source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison/DEGs/ct_DEGs_func.R")

main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs/outputs/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/"

sub_disease <- list.dirs(main_DISCO, full.names = F, recursive = F)
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

disco <- ImportDataset(main_DISCO, sub_disease)
UCSC <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes")


common_annotation <- c("CXCL14 IN" = "Interneurons",
              "EC" = "Endothelial cells",
              "fibrous astrocyte"  = "Astrocytes",
              "L2_3 EN" = "Excitatory neurons", 
              "L4 EN" = "Excitatory neurons",
              "L5 EN" = "Excitatory neurons",
              "L5_6 EN" = "Excitatory neurons",
              "L5b EN" = "Excitatory neurons",
              "L6 EN" = "Excitatory neurons",                
              "microglia" = "Microglia", 
              "oligodendrocyte" =  "Oligodendrocytes",      
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
              "Oligodendrocytes" = "Oligodendrocyte",
              "OPCs" = "OPCs",            
              "Unknown" = "Unknown",           
             "Vascular cells" = "Vascular cells",     
             "Dorsal progenitors"  = "Dorsal progenitors" ,   
             "Ventral progenitors" = "Ventral progenitors")
names(common_annotation) <- tolower(names(common_annotation))

sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]))

ggplot(sexes[["OPCs"]], aes(age, gene_id, fill=presence)) +
  geom_tile(color="#D3D3D3") +
  coord_fixed() +
  #facet_wrap(~sex, nrow = 2) +
  labs(x="Developmental Ages", y="Genes", fill="Genes found") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.spacing.x=unit(0, "lines"),
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
        axis.ticks.x=element_blank(),
        axis.text.y = element_blank(),
        legend.position = "right", 
        legend.title = element_text(size=12, face="bold", colour = "black"))