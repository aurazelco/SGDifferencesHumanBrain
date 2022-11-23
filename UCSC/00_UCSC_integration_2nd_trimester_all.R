library(Seurat)
library(patchwork)
library(metap)
library(stringr)

main <- "/Home/ii/auraz/data/UCSC/outputs/Eze_Nowakowski_Velmeshev_2nd_trimester_integrated"

dir.create(main, recursive = T, showWarnings = F)

# Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html

input_rds_path <-  "/Home/ii/auraz/data/UCSC/Seurat_UCSC/integrated"
input_rds_files <- c("/Home/ii/auraz/data/UCSC/Seurat_UCSC/integrated/Eze_Nowakowski_integrated.rds",
                     "/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_2nd_trimester.rds")
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- input_rds_files
names(input_rds) <- str_remove_all(names(input_rds), ".rds")
names(input_rds) <- str_extract(names(input_rds), "\\w+$")

# normalize and identify variable features for each dataset independently
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = input_rds)
common.anchors <- FindIntegrationAnchors(object.list = input_rds, anchor.features = features)

# this command creates an 'integrated' data assay
trim_2nd <- IntegrateData(anchorset = common.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(velm) <- "integrated"

# Run the standard workflow for visualization and clustering
trim_2nd <- ScaleData(trim_2nd, verbose = FALSE)
trim_2nd <- RunPCA(trim_2nd, npcs = 30, verbose = FALSE)
trim_2nd <- RunUMAP(trim_2nd, reduction = "pca", dims = 1:30)
trim_2nd <- FindNeighbors(trim_2nd, reduction = "pca", dims = 1:30)
trim_2nd <- FindClusters(trim_2nd, resolution = 0.5)

trim_2nd@project.name <- "trim_2nd_all"

# Visualization
pdf(paste0(main, trim_2nd@project.name, "_age.pdf"))
DimPlot(trim_2nd, reduction = "umap", group.by = "age")
dev.off()

pdf(paste0(main, trim_2nd@project.name, "_cluster_final.pdf"))
DimPlot(trim_2nd, reduction = "umap", group.by = "cluster_final")
dev.off()

saveRDS(trim_2nd, paste0(input_rds_path, trim_2nd@project.name, ".rds"))
