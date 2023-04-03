library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggpubr)

out_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/suppl_files/"
plot_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Extra_figures/"


####################################################################################################
#
# NUM OF CELLS
#
####################################################################################################


disco <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/num_proj_sex_ct.csv")

velm <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_sex_ct_per_age.csv")
#eze_nowa <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Eze_Nowakowski_num_cells.csv")

disco$X <- NULL
disco$disease <- str_replace_all(disco$disease, "Normal", "Healthy")
velm$X <- NULL
#eze_nowa$X <- NULL
#eze_nowa <- cbind(rep("Eze_Nowakowski_2nd_trimester", nrow(eze_nowa)), eze_nowa)

#colnames(eze_nowa) <- c("proj", "sex", "ct", "count")
colnames(velm) <- c("proj", "sex", "ct", "count")
#eze_nowa$disease <- rep("Normal", nrow(eze_nowa))
velm$disease <- rep("Healthy", nrow(velm))

#eze_nowa <- eze_nowa %>% relocate(disease, .after = sex)
velm <- velm %>% relocate(disease, .after = sex)

velm$proj <- paste("Velmeshev_2022", velm$proj, sep = "_")
disco$proj <- paste("DISCO", disco$proj, sep = "_")


#all_ds_ct <- rbind(eze_nowa, velm, disco)
all_ds_ct <- rbind(velm, disco)
all_ds_ct$sex <- str_replace_all(all_ds_ct$sex, c("Female"="F", "Male"="M"))

####################################################################################################
#
# NUM OF DEGS
#
####################################################################################################

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison_adjust_pval/DEGs/Compare_DEGs_func.R")

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
sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)
num_deg <- NumDEGsAcrossConditions(sexes, groups_order)

####################################################################################################
#
# CORRELATION
#
####################################################################################################

num_deg$id <- paste(num_deg$condition, num_deg$sex, num_deg$ct, sep="--")

all_degs <- num_deg[, c("id", "count_degs")]

all_ds_ct$condition <- paste0(all_ds_ct$disease, all_ds_ct$proj)
all_ds_ct$condition <- str_remove_all(all_ds_ct$condition, "DISCO")
all_ds_ct$condition <- str_replace_all(all_ds_ct$condition, c("HealthyVelmeshev"="Velmeshev", "-"="_", " "="_"))
all_ds_ct$id <- paste(all_ds_ct$condition, all_ds_ct$sex, all_ds_ct$ct, sep="--")

all_cells <- all_ds_ct[, c("id", "count")]

corr_df <- merge(all_cells, all_degs, by="id")
colnames(corr_df) <- c("id", "cells", "degs")

corr_df <- separate(corr_df, id, into = c("dataset", "sex", "ct"), sep="--")

pdf(paste0(plot_path, "corr_number_cells_vs_degs.pdf"), width = 8, height = 6)
print(
  ggplot(corr_df[which(corr_df$cells>=100), ], aes(log10(cells), degs, color=sex)) +
    geom_point(stat="identity") +
    facet_wrap( ~ sex, scales = "free") +
    labs(x="Log10 of cell number", y="Number of DEGs", color="Sex") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size=12, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black"),
          axis.text.y = element_text(size=8, colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "right", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()
