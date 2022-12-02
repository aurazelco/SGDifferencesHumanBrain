library(Seurat)
library(SeuratObject)

disco_filt <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/brainV1.0_all_FM_filt.rds")
disco_path <- "/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/"

meta_disco <-disco_filt@meta.data
mat_disco <- as.matrix(GetAssayData(disco_filt[["RNA"]], slot="data"))

rm(disco_filt)
#colnames(mat_disco) <- rownames(meta_disco)

disco <- CreateSeuratObject(counts = mat_disco, project = "DISCO", meta.data=meta_disco, assay = "RNA")

rm(meta_disco, mat_disco, genes)

disco[["percent.mt"]] <- PercentageFeatureSet(disco, pattern = "^MT-")
disco <- subset(disco, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
disco <- NormalizeData(disco)
disco <- FindVariableFeatures(disco, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(disco)
disco <- ScaleData(disco, features = all.genes)
disco <- RunPCA(disco, features = VariableFeatures(object = disco))
VizDimLoadings(disco, dims = 1:2, reduction = "pca")
disco <- FindNeighbors(disco, dims = 1:15)
disco <- FindClusters(disco, resolution = 0.5)
disco <- RunUMAP(disco, dims = 1:15)

