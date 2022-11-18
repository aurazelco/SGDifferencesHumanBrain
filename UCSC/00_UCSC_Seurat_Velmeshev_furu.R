library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)

out_path <- "/Home/ii/auraz/data/UCSC"

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
# 2ND TRIMESTER
#
####################################################################################################

meta <- read.csv(paste0(out_path, "/UCSC_downloads/new_meta_Velmeshev_2022.csv"), header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell

meta_2nd_trim <- subset(meta, age=="2nd trimester")

mat_2nd_trim_1 <- fread(paste0(out_path, "/UCSC_downloads/2nd_trimester_1.tsv.gz"))
mat_2nd_trim_2 <- fread(paste0(out_path, "/UCSC_downloads/2nd_trimester_2.tsv.gz"))
mat_2nd_trim_2[,1] <- NULL
colnames(mat_2nd_trim_1) <- c("gene", (seq(1,(ncol(mat_2nd_trim_1)-1))))
colnames(mat_2nd_trim_2) <- c("gene", (seq(ncol(mat_2nd_trim_1), (ncol(mat_2nd_trim_1) + ncol(mat_2nd_trim_2) - 2))))
mat_2nd_trim <- merge(mat_2nd_trim_1, mat_2nd_trim_2, by="gene")
rm(mat_2nd_trim_1, mat_2nd_trim_2)

genes <- mat_2nd_trim[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_2nd_trim <- data.frame(mat_2nd_trim[,-1], row.names=genes)
colnames(mat_2nd_trim) <- rownames(meta_2nd_trim)

velm_2nd_trim <- CreateSeuratObject(counts = mat_2nd_trim, project = "Velmeshev_2022_2nd_trimester", meta.data=meta_2nd_trim, assay = "RNA")
rm(meta_2nd_trim, mat_2nd_trim)

#saveRDS(velm_2nd_trim, paste0(out_path, "/Seurat_UCSC/", velm_2nd_trim@project.name, ".rds"))

#velm_2nd_trim <- readRDS("Seurat_UCSC/Velmeshev_2022_2nd_trimester.rds")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_2nd_trim[["percent.mt"]] <- PercentageFeatureSet(velm_2nd_trim, pattern = "^MT-")
# Visualize QC metrics as a violin plot
#VlnPlot(velm_2nd_trim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(velm_2nd_trim, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(velm_2nd_trim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
velm_2nd_trim <- subset(velm_2nd_trim, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_2nd_trim <- NormalizeData(velm_2nd_trim)
velm_2nd_trim <- FindVariableFeatures(velm_2nd_trim, selection.method = "vst", nfeatures = 2000)
#saveRDS(velm_2nd_trim, paste0("Seurat_UCSC/", velm_2nd_trim@project.name, ".rds"))
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(velm_2nd_trim), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(velm_2nd_trim)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(velm_2nd_trim)
velm_2nd_trim <- ScaleData(velm_2nd_trim, features = all.genes)
velm_2nd_trim <- RunPCA(velm_2nd_trim, features = VariableFeatures(object = velm_2nd_trim))
# Examine and visualize PCA results a few different ways
#print(velm_2nd_trim[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_2nd_trim, dims = 1:2, reduction = "pca")
#DimPlot(velm_2nd_trim, reduction = "pca") 
#DimHeatmap(velm_2nd_trim, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#velm_2nd_trim <- JackStraw(velm_2nd_trim, num.replicate = 100)
#velm_2nd_trim <- ScoreJackStraw(velm_2nd_trim, dims = 1:20)
#JackStrawPlot(velm_2nd_trim, dims = 1:15)
#ElbowPlot(velm_2nd_trim)
# based on rprevious plots, decide the number of dimensions
velm_2nd_trim <- FindNeighbors(velm_2nd_trim, dims = 1:15)
velm_2nd_trim <- FindClusters(velm_2nd_trim, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_2nd_trim), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_2nd_trim <- RunUMAP(velm_2nd_trim, dims = 1:15)
#DimPlot(velm_2nd_trim, reduction = "umap")

velm_2nd_trim@meta.data$sex_age <- paste(velm_2nd_trim@meta.data$sex, velm_2nd_trim@meta.data$age, sep="_")
velm_2nd_trim@meta.data$proj <-  rep(velm_2nd_trim@project.name, nrow(velm_2nd_trim@meta.data))
velm_2nd_trim@meta.data$id_sex_age <- paste(velm_2nd_trim@meta.data$samples, velm_2nd_trim@meta.data$sex, velm_2nd_trim@meta.data$age, sep="_")

#pdf(paste0(out_path, "/outputs/", velm_2nd_trim@project.name,  "_XIST.pdf"))
#print(VlnPlot(velm_2nd_trim, features = "XIST", group.by = "id_sex_age") + NoLegend())
#dev.off()
# XIST and other Ygenes from 00_filterin_disco not present -> we have to trust the metadata

saveRDS(velm_2nd_trim, paste0(out_path, "/Seurat_UCSC/", velm_2nd_trim@project.name, ".rds"))

velm_num_cells <- as.data.frame(table(velm_2nd_trim$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_2nd_trim@project.name, nrow(velm_num_cells)), velm_num_cells)
write.csv(velm_num_cells, paste0(out_path, "/outputs/", velm_2nd_trim@project.name, "_num_cells.csv"))

################

velm_2nd_trim@meta.data$cluster_final <- rep("no_data", nrow(velm_2nd_trim@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_2nd_trim@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_2nd_trim@meta.data[which(velm_2nd_trim@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}

pdf(paste0(out_path, "/outputs/", velm_2nd_trim@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_2nd_trim, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_2nd_trim, paste0(out_path, "/Seurat_UCSC/", velm_2nd_trim@project.name, ".rds"))


#velm_2nd_trim <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev_2022_2nd_trimester.rds")


################

velm_2nd_trim <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev_2022_2nd_trimester.rds")

velm_2nd_trim@meta.data$sex_ct <- paste(velm_2nd_trim@meta.data$sex, velm_2nd_trim@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_2nd_trim$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_2nd_trim@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

write.csv(velm_num_sex_ct, paste0(out_path, "/outputs/", velm_2nd_trim@project.name, "_num_sex_ct_per_age.csv"))


####################################################################################################
#
# 10-20 YEARS
#
####################################################################################################

meta_10_20_years <- subset(meta, age=="10-20 years")

mat_10_20_years <- fread(paste0(out_path, "/UCSC_downloads/10_20_years.tsv.gz"))
genes <- mat_10_20_years[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_10_20_years <- data.frame(mat_10_20_years[,-1], row.names=genes)
colnames(mat_10_20_years) <- rownames(meta_10_20_years)

velm_10_20_years <- CreateSeuratObject(counts = mat_10_20_years, project = "Velmeshev_2022_10_20_years", meta.data=meta_10_20_years, assay = "RNA")

rm(meta_10_20_years, mat_10_20_years, genes)

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_10_20_years[["percent.mt"]] <- PercentageFeatureSet(velm_10_20_years, pattern = "^MT-")
# Visualize QC metrics as a violin plot
#VlnPlot(velm_10_20_years, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(velm_10_20_years, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(velm_10_20_years, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
velm_10_20_years <- subset(velm_10_20_years, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_10_20_years <- NormalizeData(velm_10_20_years)
velm_10_20_years <- FindVariableFeatures(velm_10_20_years, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
#top10 <- head(VariableFeatures(velm_10_20_years), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(velm_10_20_years)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(velm_10_20_years)
velm_10_20_years <- ScaleData(velm_10_20_years, features = all.genes)
velm_10_20_years <- RunPCA(velm_10_20_years, features = VariableFeatures(object = velm_10_20_years))
# Examine and visualize PCA results a few different ways
#print(velm_10_20_years[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_10_20_years, dims = 1:2, reduction = "pca")
#DimPlot(velm_10_20_years, reduction = "pca") 
#DimHeatmap(velm_10_20_years, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
#velm_10_20_years <- JackStraw(velm_10_20_years, num.replicate = 100)
#velm_10_20_years <- ScoreJackStraw(velm_10_20_years, dims = 1:20)
#JackStrawPlot(velm_10_20_years, dims = 1:15)
#ElbowPlot(velm_10_20_years)
# based on rprevious plots, decide the number of dimensions
velm_10_20_years <- FindNeighbors(velm_10_20_years, dims = 1:15)
velm_10_20_years <- FindClusters(velm_10_20_years, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_10_20_years), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_10_20_years <- RunUMAP(velm_10_20_years, dims = 1:15)
#DimPlot(velm_10_20_years, reduction = "umap")

velm_10_20_years@meta.data$sex_age <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$age, sep="_")
velm_10_20_years@meta.data$proj <-  rep(velm_10_20_years@project.name, nrow(velm_10_20_years@meta.data))
velm_10_20_years@meta.data$id_sex_age <- paste(velm_10_20_years@meta.data$samples, velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$age, sep="_")

pdf(paste0(out_path, "/outputs/", velm_10_20_years@project.name,  "_XIST.pdf"))
print(VlnPlot(velm_10_20_years, features = "XIST", group.by = "id_sex_age") + NoLegend())
dev.off()

#saveRDS(velm_10_20_years, paste0(out_path, "/Seurat_UCSC/", velm_10_20_years@project.name, ".rds"))

velm_num_cells <- as.data.frame(table(velm_10_20_years$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_10_20_years@project.name, nrow(velm_num_cells)), velm_num_cells)
write.csv(velm_num_cells, paste0(out_path, "/outputs/", velm_10_20_years@project.name, "_num_cells.csv"))


################

velm_10_20_years@meta.data$cluster_final <- rep("no_data", nrow(velm_10_20_years@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_10_20_years@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_10_20_years@meta.data[which(velm_10_20_years@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}

pdf(paste0(out_path, "/outputs/", velm_10_20_years@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_10_20_years, reduction = "umap", group.by = "cluster_final"))
dev.off()


saveRDS(velm_10_20_years, paste0(out_path, "/Seurat_UCSC/", velm_10_20_years@project.name, ".rds"))


#velm_10_20_years <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev_2022_10_20_years.rds")

################


velm_10_20_years <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev_2022_10_20_years.rds")

velm_10_20_years@meta.data$sex_ct <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_10_20_years$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_10_20_years@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

write.csv(velm_num_sex_ct, paste0(out_path, "/outputs/", velm_10_20_years@project.name, "_num_sex_ct_per_age.csv"))
