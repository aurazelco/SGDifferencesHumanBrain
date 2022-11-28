library(Seurat)
library(SeuratObject)
library(stringr)

main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/"

rds_path <-  "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC"

rds_files <- list.files(rds_path, "Velmeshev_2022", full.names = T)

names(rds_files) <- str_remove_all(rds_files, paste0(rds_path, "/Velmeshev_2022_"))
names(rds_files) <- str_remove_all(names(rds_files), ".rds")

##### 

####################################################################################################
#
# 3RD TRIMESTER
#
####################################################################################################

velm_rds <- readRDS(rds_files[["3rd_trimester"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))
rm(velm_rds)


####################################################################################################
#
# 0-1 YEARS
#
####################################################################################################

velm_rds <- readRDS(rds_files[["0_1_years"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))
rm(velm_rds)

####################################################################################################
#
# 1-2- YEARS
#
####################################################################################################

velm_rds <- readRDS(rds_files[["1_2_years"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))

VlnPlot(velm_rds, features = "XIST", group.by = "id_sex_age") + NoLegend()

# one sample is actually F and not M -> chnage that into metadata as well

velm_rds@meta.data[which(velm_rds@meta.data$samples=="1-1547-BA24"), "sex"] <- "Female"
velm_rds@meta.data$sex_age <- paste(velm_rds@meta.data$sex, velm_rds@meta.data$age, sep="_")
velm_rds@meta.data$proj <-  rep(velm_rds@project.name, nrow(velm_rds@meta.data))
velm_rds@meta.data$id_sex_age <- paste(velm_rds@meta.data$samples, velm_rds@meta.data$sex, velm_rds@meta.data$age, sep="_")


rm(velm_rds)
####################################################################################################
#
# 2-4 YEARS
#
####################################################################################################

velm_rds <- readRDS(rds_files[["2_4_years"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))
rm(velm_rds)

####################################################################################################
#
# ADULT
#
####################################################################################################

velm_rds <- readRDS(rds_files[["Adult"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))
rm(velm_rds)