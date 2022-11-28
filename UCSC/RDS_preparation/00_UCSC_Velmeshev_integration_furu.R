library(Seurat)
library(patchwork)
library(stringr)

main <- "/Home/ii/auraz/data/UCSC/outputs/integrated"


# Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html
input_rds_path <-  "/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev"
input_rds_files <- list.files(path = input_rds_path, pattern = ".rds", full.names = T)
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- list.files(path = input_rds_path, pattern = ".rds", full.names = F)
names(input_rds) <- str_remove_all(names(input_rds), ".rds")
names(input_rds) <- str_remove_all(names(input_rds), "Velmeshev_2022_")

# normalize and identify variable features for each dataset independently
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = input_rds)
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# this command creates an 'integrated' data assay
anchors <- FindIntegrationAnchors(object.list = input_rds, reduction = "rpca", dims = 1:50)
velm <- IntegrateData(anchorset = anchors, dims = 1:50)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(velm) <- "integrated"

# Run the standard workflow for visualization and clustering
velm <- ScaleData(velm, verbose = FALSE)
velm <- RunPCA(velm, npcs = 30, verbose = FALSE)
velm <- RunUMAP(velm, reduction = "pca", dims = 1:30)
velm <- FindNeighbors(velm, reduction = "pca", dims = 1:30)
velm <- FindClusters(velm, resolution = 0.5)

velm@project.name <- "Velmeshev_all"

# Visualization
pdf(paste0(main, velm@project.name, "_age.pdf"))
DimPlot(velm, reduction = "umap", group.by = "age")
dev.off()

pdf(paste0(main, velm@project.name, "_cluster_final.pdf"))
DimPlot(velm, reduction = "umap", group.by = "cluster_final")
dev.off()

saveRDS(velm, paste0(input_rds_path, velm@project.name, ".rds"))
