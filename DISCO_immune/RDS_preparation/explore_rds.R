library(Seurat)
library(SeuratObject)

rds_path <- "/Home/ii/auraz/data/DISCO_immune/"

## BONE MARROW

bm_folder <- paste0(rds_path, "Bone_marrow/")
bm_Zenodo <- readRDS(paste0(rds_path, "disco_bone_marrow_v2.0.rds"))

bm_exprMtx <- GetAssayData(bm_Zenodo[["RNA"]], slot="counts")

bm <- CreateSeuratObject(counts = bm_exprMtx, project = "DISCO_bone_marrow", meta.data=bm_Zenodo@meta.data, assay = "RNA")

rm(bm_Zenodo, bm_exprMtx)

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
bm[["percent.mt"]] <- PercentageFeatureSet(bm, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(bm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(bm, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bm, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
bm <- subset(bm, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
bm <- NormalizeData(bm)
bm <- FindVariableFeatures(bm, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(bm), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(bm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(bm)
bm <- ScaleData(bm, features = all.genes)
bm <- RunPCA(bm, features = VariableFeatures(object = bm))
# Examine and visualize PCA results a few different ways
print(bm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(bm, dims = 1:2, reduction = "pca")
DimPlot(bm, reduction = "pca") 
DimHeatmap(bm, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
bm <- JackStraw(bm, num.replicate = 100)
bm <- ScoreJackStraw(bm, dims = 1:20)
JackStrawPlot(bm, dims = 1:15)
ElbowPlot(bm)
# based on rprevious plots, decide the number of dimensions
bm <- FindNeighbors(bm, dims = 1:12)
bm <- FindClusters(bm, resolution = 0.5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
bm <- RunUMAP(bm, dims = 1:12)


## THYMUS

thy_folder <- paste0(rds_path, "Thymus/")
thy_Zenodo <- readRDS(paste0(rds_path, "disco_thymus_v1.0.rds"))

thy_exprMtx <- GetAssayData(thy_Zenodo[["RNA"]], slot="counts")

thy <- CreateSeuratObject(counts = thy_exprMtx, project = "DISCO_thymus", meta.data=thy_Zenodo@meta.data, assay = "RNA")

rm(thy_Zenodo, thy_exprMtx)

# Seurat tutorial https://satijalab.org/seurat/articles/pthyc3k_tutorial.html
thy[["percent.mt"]] <- PercentageFeatureSet(thy, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(thy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(thy, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(thy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
thy <- subset(thy, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
thy <- NormalizeData(thy)
thy <- FindVariableFeatures(thy, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(thy), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(thy)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(thy)
thy <- ScaleData(thy, features = all.genes)
thy <- RunPCA(thy, features = VariableFeatures(object = thy))
# Examine and visualize PCA results a few different ways
print(thy[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(thy, dims = 1:2, reduction = "pca")
DimPlot(thy, reduction = "pca") 
DimHeatmap(thy, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
thy <- JackStraw(thy, num.replicate = 100)
thy <- ScoreJackStraw(thy, dims = 1:20)
JackStrawPlot(thy, dims = 1:15)
ElbowPlot(thy)
# based on rprevious plots, decide the number of dimensions
thy <- FindNeighbors(thy, dims = 1:12)
thy <- FindClusters(thy, resolution = 0.5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
thy <- RunUMAP(thy, dims = 1:12)

