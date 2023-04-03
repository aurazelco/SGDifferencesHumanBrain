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

VlnPlot(velm_rds, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(main, "TMSB4_RNA_expression_",velm_rds@project.name, ".pdf"), width = 10, height = 8)

rm(velm_rds)


####################################################################################################
#
# 0-1 YEARS
#
####################################################################################################

velm_rds <- readRDS(rds_files[["0_1_years"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))

VlnPlot(velm_rds, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(main, "TMSB4_RNA_expression_",velm_rds@project.name, ".pdf"), width = 10, height = 8)

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

# one sample is actually F and not M -> change that into metadata as well

velm_rds@meta.data[which(velm_rds@meta.data$samples=="1-1547-BA24"), "sex"] <- "Female"
velm_rds@meta.data$sex_age <- paste(velm_rds@meta.data$sex, velm_rds@meta.data$age, sep="_")
velm_rds@meta.data$proj <-  rep(velm_rds@project.name, nrow(velm_rds@meta.data))
velm_rds@meta.data$id_sex_age <- paste(velm_rds@meta.data$samples, velm_rds@meta.data$sex, velm_rds@meta.data$age, sep="_")

VlnPlot(velm_rds, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(main, "TMSB4_RNA_expression_",velm_rds@project.name, ".pdf"), width = 10, height = 8)


rm(velm_rds)

####################################################################################################
#
# 2-4 YEARS
#
####################################################################################################

velm_rds <- readRDS(rds_files[["2_4_years"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))

VlnPlot(velm_rds, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(main, "TMSB4_RNA_expression_",velm_rds@project.name, ".pdf"), width = 10, height = 8)

rm(velm_rds)

####################################################################################################
#
# ADULT
#
####################################################################################################

velm_rds <- readRDS(rds_files[["Adult"]])
DimPlot(velm_rds, reduction="umap")
ggsave(paste0(main, velm_rds@project.name, "_UMAP_cluster.pdf"))

VlnPlot(velm_rds, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(main, "TMSB4_RNA_expression_",velm_rds@project.name, ".pdf"), width = 10, height = 8)

rm(velm_rds)


####################################################################################################
#
# 2ND TRIMESTER
#
####################################################################################################

out_path <- "/Home/ii/auraz/data/UCSC"
rds_path <- paste0(out_path, "/Seurat_UCSC/Velmeshev/")
output_2nd_trim <- paste0(out_path, "/outputs/Velmeshev_2nd_trimester/")

velm_2nd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_2nd_trimester.rds"))

VlnPlot(velm_2nd_trim, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(output_2nd_trim, "TMSB4_RNA_expression_",velm_2nd_trim@project.name, ".pdf"), width = 10, height = 8)

rm(velm_2nd_trim)



####################################################################################################
#
# 10-20 YEARS
#
####################################################################################################
output_10_20_yo <- paste0(out_path, "/outputs/Velmeshev_10_20_years")

velm_10_20_years <- readRDS(paste0(rds_path, "Velmeshev_2022_10_20_years.rds"))

VlnPlot(velm_10_20_years, features = c("TMSB4X", "TMSB4Y"), group.by = "sex" ) 
ggsave(paste0(output_10_20_yo, "TMSB4_RNA_expression_", velm_10_20_years@project.name, ".pdf"), width = 10, height = 8)

rm(velm_10_20_years)
