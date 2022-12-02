library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)

# Clusters annotation using Fig 1B and Fig S1a from paper

ann_clusters <- list("Astrocytes" = c(3),
                     "Microglia" = c(25),
                     "OPCs" = c(9),
                     "Oligodendrocytes" = c(11),
                     "Excitatory neurons" = c(0,13,4,7,21,12,19,24,6,28,23,10,22,26,15, 30),
                     "Dorsal progenitors" = c(18),
                     "Ventral progenitors" = c(5),
                     "Interneurons" = c(17, 31,32,8,16, 14),
                     "Vascular cells" = c(20),
                     "Debris" = c(1,2),
                     "Unknown" = c(29, 27)
)

Reduce(intersect, ann_clusters)

cts <- vector()
og_clusters <- vector()

for (ct in names(ann_clusters)) {
  cts <- c(cts, rep(ct, length(ann_clusters[[ct]])))
  og_clusters <- c(og_clusters, ann_clusters[[ct]])
}
ann_df <- data.frame(cts, og_clusters)

rm(ann_clusters)

####################################################################################################
#
# ALL AGES
#
####################################################################################################

out_path <- "/Home/ii/auraz/data/UCSC/outputs/Velmeshev_integrated/"
rds_path <- "/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/"
dir.create(out_path, recursive=T, showWarnings = F)
dir.create(rds_path, recursive=T, showWarnings = F)

meta <- read.csv("/Home/ii/auraz/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", 
                 header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell

mat_Velm_all <- fread("/Home/ii/auraz/data/UCSC/UCSC_downloads/exprMtx_filt_Velmeshev_2022.tsv.gz")
genes <- mat_Velm_all[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_Velm_all <- data.frame(mat_Velm_all[,-1], row.names=genes)
colnames(mat_Velm_all) <- rownames(meta)

Velm_all <- CreateSeuratObject(counts = mat_Velm_all, project = "Velmeshev_2022", meta.data=meta, assay = "RNA")
rm(meta, mat_Velm_all)

#saveRDS(Velm_all, paste0(out_path, "/Seurat_UCSC/", Velm_all@project.name, ".rds"))

#Velm_all <- readRDS("Seurat_UCSC/Velmeshev_2022_2nd_trimester.rds")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
Velm_all[["percent.mt"]] <- PercentageFeatureSet(Velm_all, pattern = "^MT-")
# Visualize QC metrics as a violin plot
#VlnPlot(Velm_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(Velm_all, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(Velm_all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
Velm_all <- subset(Velm_all, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
Velm_all <- NormalizeData(Velm_all)
Velm_all <- FindVariableFeatures(Velm_all, selection.method = "vst", nfeatures = 2000)
#saveRDS(Velm_all, paste0("Seurat_UCSC/", Velm_all@project.name, ".rds"))
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(Velm_all), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(Velm_all)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(Velm_all)
Velm_all <- ScaleData(Velm_all, features = all.genes)
Velm_all <- RunPCA(Velm_all, features = VariableFeatures(object = Velm_all))
# Examine and visualize PCA results a few different ways
#print(Velm_all[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Velm_all, dims = 1:2, reduction = "pca")
#DimPlot(Velm_all, reduction = "pca") 
#DimHeatmap(Velm_all, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#Velm_all <- JackStraw(Velm_all, num.replicate = 100)
#Velm_all <- ScoreJackStraw(Velm_all, dims = 1:20)
#JackStrawPlot(Velm_all, dims = 1:15)
#ElbowPlot(Velm_all)
# based on rprevious plots, decide the number of dimensions
Velm_all <- FindNeighbors(Velm_all, dims = 1:15)
Velm_all <- FindClusters(Velm_all, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(Velm_all), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
Velm_all <- RunUMAP(Velm_all, dims = 1:15)
#DimPlot(Velm_all, reduction = "umap")

print(colnames(Velm_all@meta.data))

Velm_all@meta.data$sex_age <- paste(Velm_all@meta.data$sex, Velm_all@meta.data$age, sep="_")

#saveRDS(Velm_all, paste0(rds_path, Velm_all@project.name, ".rds"))

################

Velm_all@meta.data$cluster_final <- rep("no_data", nrow(Velm_all@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% Velm_all@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    Velm_all@meta.data[which(Velm_all@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}

pdf(paste0(out_path, Velm_all@project.name,  "_cluster_final.pdf"))
print(DimPlot(Velm_all, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(Velm_all, paste0(rds_path, Velm_all@project.name, ".rds"))


################

pdf(paste0(out_path, Velm_all@project.name,  "_age.pdf"))
print(DimPlot(Velm_all, reduction = "umap", group.by = "age"))
dev.off()

pdf(paste0(out_path, Velm_all@project.name,  "_sex.pdf"))
print(DimPlot(Velm_all, reduction = "umap", group.by = "sex"))
dev.off()

pdf(paste0(out_path, Velm_all@project.name,  "_sex_age.pdf"))
print(DimPlot(Velm_all, reduction = "umap", group.by = "sex_age"))
dev.off()
