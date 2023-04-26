# Author: Aura Zelco
# Brief description:
  # This script is used to prepare supplemental files for the paper
  # Brief procedure: 
  # 1. generates plots for XIST expression in samples
  # 2. generates plots for TMSB4X and TMSB4Y expression in sexes

# OBS: this script requires visual inspection and manually-entered information, therefore should be opened in a R environment/IDE (e.g. RStudio). 

# Imports necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)

# Defines the main directories
disco_path <- "DISCOv1.0/"


##################################### 1. XIST expression

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

disco_filt@meta.data$id_sex <- paste(disco_filt@meta.data$sample_id, disco_filt@meta.data$gender, sep = "_")

VlnPlot(disco_filt, features = "XIST", group.by = "id_sex" ) + NoLegend()
ggsave(paste0(disco_path, "XIST_RNA_expression.pdf"))

##################################### 2. TMSB4X and TMSB4Y expression

VlnPlot(disco_filt, features = c("TMSB4X", "TMSB4Y"), group.by = "gender" ) 
ggsave(paste0(disco_path, "TMSB4_RNA_expression.pdf"), width = 10, height = 8)
