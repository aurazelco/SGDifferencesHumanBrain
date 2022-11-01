require(Seurat)
require(data.table)

main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/NEMO_fetal/outputs/"

# from this tutorial UCSC CellBrowser https://cellbrowser.readthedocs.io/en/master/load.html
mat2 <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/NEMO_fetal/exprMatrix.tsv.gz")

mat <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/NEMO_fetal/counts_exprMatrix.tsv.gz")
meta <- read.table("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/NEMO_fetal/meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
genes <- mat[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[,-1], row.names=genes)
nemo <- CreateSeuratObject(counts = mat, project = "NEMOfetal", meta.data=meta, assay = "RNA")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
nemo[["percent.mt"]] <- PercentageFeatureSet(nemo, pattern = "^MT-")
# Visualize QC metrics as a violin plot
#VlnPlot(nemo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(nemo, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(nemo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
#nemo <- subset(nemo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
nemo <- subset(nemo, subset = percent.mt < 5)
nemo <- NormalizeData(nemo)
nemo <- FindVariableFeatures(nemo, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nemo), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(nemo)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(nemo)
nemo <- ScaleData(nemo, features = all.genes)
nemo <- RunPCA(nemo, features = VariableFeatures(object = nemo))
# Examine and visualize PCA results a few different ways
print(nemo[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(nemo, dims = 1:2, reduction = "pca")
DimPlot(nemo, reduction = "pca") 
DimHeatmap(nemo, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
nemo <- JackStraw(nemo, num.replicate = 100)
nemo <- ScoreJackStraw(nemo, dims = 1:20)
p2 <- JackStrawPlot(nemo, dims = 1:15)
p3 <- ElbowPlot(nemo)
plot<- ggarrange(p2, p3, labels = c("JackStrawPlot", "ElbowPlot"), ncol = 2)
print(
  annotate_figure(plot, top = text_grob(id, 
                                        color = "black", face = "bold"))
)
id_seurat <- FindNeighbors(id_seurat, dims = 1:cluster_num)
id_seurat <- FindClusters(id_seurat, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(id_seurat), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
id_seurat <- RunUMAP(id_seurat, dims = 1:cluster_num)