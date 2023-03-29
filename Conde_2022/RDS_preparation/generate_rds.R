library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)


files_path <- "/Home/ii/auraz/data/Conde_2022/"

meta <- read.csv(paste0(files_path, "original_files/meta.tsv"), sep="\t")

mat <- fread(paste0(files_path, "original_files/exprMatrix.tsv.gz"))
genes <- mat[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[,-1], row.names=genes)
rownames(meta) <- meta$cellId
colnames(mat) <- rownames(meta)

conde <- CreateSeuratObject(counts = mat, project = "Conde_2022", meta.data=meta, assay = "RNA")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
conde[["percent.mt"]] <- PercentageFeatureSet(conde, pattern = "^MT-")
# Visualize QC metrics as a violin plot
#VlnPlot(conde, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
#plot1 <- FeatureScatter(conde, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(conde, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2
conde <- subset(conde, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
saveRDS(conde, paste0(files_path, "conde_2022.rds"))

conde <- NormalizeData(conde)
conde <- FindVariableFeatures(conde, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(conde), 10)
# plot variable features with and without labels
#plot1 <- VariableFeaturePlot(conde)
#plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#plot1 + plot2
all.genes <- rownames(conde)
conde <- ScaleData(conde, features = all.genes)
conde <- RunPCA(conde, features = VariableFeatures(object = conde))
# Examine and visualize PCA results a few different ways
#print(conde[["pca"]], dims = 1:5, nfeatures = 5)
#VizDimLoadings(conde, dims = 1:2, reduction = "pca")
#DimPlot(conde, reduction = "pca") 
#DimHeatmap(conde, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
conde <- JackStraw(conde, num.replicate = 100)
conde <- ScoreJackStraw(conde, dims = 1:20)
#JackStrawPlot(conde, dims = 1:15)
#ElbowPlot(conde)
# based on rprevious plots, decide the number of dimensions
conde <- FindNeighbors(conde, dims = 1:12)
conde <- FindClusters(conde, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(conde), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
conde <- RunUMAP(conde, dims = 1:12)

conde@meta.data$Sex <- str_replace_all(conde@meta.data$Sex, c("Female"="F", "Male"="M"))
saveRDS(conde, paste0(files_path, "conde_2022.rds"))

out_path <- paste0(files_path, "outputs/")
dir.create(out_path, showWarnings = F, recursive = T)

##### Saves num cells per cell type, organ and sex

conde@meta.data$donor_sex_organ_ct <- paste(conde@meta.data$Donor, conde@meta.data$Sex, conde@meta.data$Organ, conde@meta.data$Manually_curated_celltype, sep = "--")
cell_counts <- as.data.frame(table(conde@meta.data$donor_sex_organ_ct))
cell_counts <- separate(cell_counts, Var1, into = c("donor", "sex", "organ", "ct"), sep = "--")

write.csv(cell_counts, paste0(out_path, "num_cells.csv"))

pdf(paste0(out_path, "num_cells.pdf"), width = 15, height = 15)
print(
  ggplot(cell_counts, aes(ct, log10(Freq), fill=sex)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept = log10(100), linetype="dashed") +
    facet_wrap(~ organ, scales = "free") +
    labs(x="Cell types", y="Log10 cell counts", fill="Sex") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.text = element_text(size=12, colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()

pdf(paste0(out_path, "num_samples.pdf"), width = 15, height = 15)
print(
  ggplot(cell_counts, aes(donor, log10(Freq), fill=sex)) +
    geom_bar(stat="identity", position="dodge") +
    geom_hline(yintercept = log10(100), linetype="dashed") +
    facet_wrap(~ organ, scales = "free") +
    labs(x="Donors", y="Log10 cell counts", fill="Sex") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          strip.text = element_text(size=12, colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()



##### Saves UMAPs

# Sex
pdf(paste0(out_path, "umap_sex.pdf"), width = 10, height = 10)
print(
  DimPlot(conde, reduction="umap", group.by = "Sex")
)
dev.off()

pdf(paste0(out_path, "umap_donors.pdf"), width = 10, height = 10)
print(
  DimPlot(conde, reduction="umap", group.by = "Donor")
)
dev.off()

pdf(paste0(out_path, "umap_cts.pdf"), width = 10, height = 15)
print(
  DimPlot(conde, reduction="umap", group.by = "Manually_curated_celltype") + 
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "vertical")
)
dev.off()


saveRDS(conde, paste0(files_path, "conde_2022.rds"))

