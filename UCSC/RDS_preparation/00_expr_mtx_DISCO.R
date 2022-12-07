library(Seurat)
library(SeuratObject)
library(patchwork)
library(stringr)
library(dplyr)
library(data.table)

main <- "/Home/ii/auraz/data/UCSC/outputs/DISCO_Velmeshev_Adult_integrated"

dir.create(main, recursive = T, showWarnings = F)

disco <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/brainV1.0_all_FM_filt.rds")

disco_mat <- GetAssayData(disco[["RNA"]], slot="data")
disco_meta <- disco@meta.data

missing_cols <- names(which(colSums(is.na(disco_meta)) > 0))
disco_meta <- disco_meta %>% select(-missing_cols)

disco_2 <- CreateSeuratObject(counts = disco_mat, project = "DISCO", meta.data=disco_meta, assay = "RNA")

disco_projs  <- SplitObject(disco_2, split.by = "project_id")

Velmeshev <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_Adult.rds")
input_rds <- append(disco_projs, list("Velmeshev_Adult" = Velmeshev))

input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = input_rds)
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = input_rds, reduction = "rpca", dims = 1:50)
#anchors <- FindIntegrationAnchors(object.list = input_rds, reference = c(1, 2), reduction = "rpca",dims = 1:50)
adults <- IntegrateData(anchorset = anchors, dims = 1:30)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(adults) <- "integrated"

# Run the standard workflow for visualization and clustering
adults <- ScaleData(adults, verbose = T)
adults <- RunPCA(adults, npcs = 30, verbose = T)
adults <- RunUMAP(adults, reduction = "pca", dims = 1:30)
adults <- FindNeighbors(adults, reduction = "pca", dims = 1:30)
adults <- FindClusters(adults, resolution = 0.5)

adults@project.name <- "adults_all"
