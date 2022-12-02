library(Seurat)
library(SeuratObject)
library(patchwork)
library(stringr)
library(dplyr)


main <- "/Home/ii/auraz/data/UCSC/outputs/DISCO_Velmeshev_Adult_integrated"

dir.create(main, recursive = T, showWarnings = F)

input_rds_path <-  "/Home/ii/auraz/data/UCSC/Seurat_UCSC/integrated"
input_rds_files <- c("/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/brainV1.0_all_FM_filt.rds",
                     "/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_Adult.rds")
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- input_rds_files
names(input_rds) <- str_remove_all(names(input_rds), ".rds")
names(input_rds) <- str_extract(names(input_rds), "\\w+$")
names(input_rds)[1] <- "disco" 

disco_projs  <- SplitObject(input_rds[["disco"]], split.by = "project_id")

input_rds <- c(input_rds["Velmeshev_2022_Adult"], disco_projs)

input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# this command above raises some issues -> not sure how to solve it














# Next, select features for downstream integration, and run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow
features <- SelectIntegrationFeatures(object.list = input_rds)
input_rds <- lapply(X = input_rds, FUN = function(x) {
  
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# this command creates an 'integrated' data assay
anchors <- FindIntegrationAnchors(object.list = input_rds, reduction = "rpca", dims = 1:50)
#anchors <- FindIntegrationAnchors(object.list = input_rds, reference = c(1, 2), reduction = "rpca",dims = 1:50)
adults <- IntegrateData(anchorset = anchors, dims = 1:50)

DefaultAssay(adults) <- "integrated"

# Run the standard workflow for visualization and clustering
adults <- ScaleData(adults, verbose = T)
adults <- RunPCA(adults, npcs = 30, verbose = T)
adults <- RunUMAP(adults, reduction = "pca", dims = 1:30)
adults <- FindNeighbors(adults, reduction = "pca", dims = 1:30)
adults <- FindClusters(adults, resolution = 0.5)

adults@project.name <- "adults_all"

colnames(adults@meta.data)

# Visualization
pdf(paste0(main, "/", adults@project.name, "_age.pdf"))
DimPlot(adults, reduction = "umap", group.by = "age")
dev.off()

pdf(paste0(main, "/", adults@project.name, "_cluster_final.pdf"))
DimPlot(adults, reduction = "umap", group.by = "cluster_final")
dev.off()

saveRDS(adults, paste0(input_rds_path, "/", adults@project.name, ".rds"))
