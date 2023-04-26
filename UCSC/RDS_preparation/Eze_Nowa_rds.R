library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(scales)

main <- "UCSC/outputs/"


# from this tutorial UCSC CellBrowser https://cellbrowser.readthedocs.io/en/master/load.html

####################################################################################################
#
# EZE 2021
#
####################################################################################################

mat <- fread("UCSC/UCSC_downloads/exprMatrix_Eze_2021.tsv.gz")
genes <- mat[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[,-1], row.names=genes)



meta <- read.csv("UCSC/UCSC_downloads/new_meta_Eze_2021.csv", header=T, sep=",", as.is=T, row.names=1)
meta <- subset(meta, subset = Cell %in% colnames(mat))
keep <- c("Cell", "nGene.x", "nUMI.x", "Age_.Carnegie_Stage.","Individual.x","Area_As_Annotated", "Area.x", 
          "V1", "Cluster.x", "Cell.Type" ,  "Age",  "Age.Range",  "Progenitor.Cluster.Original" , "Progenitor.Combined.Cluster")

meta <- meta[, c("Cell", "nGene.x", "nUMI.x", "Age_.Carnegie_Stage.","Individual.x","Area_As_Annotated", "Area.x", 
                 "V1", "Cluster.x", "Cell.Type" ,  "Age",  "Age.Range",  "Progenitor.Cluster.Original" , "Progenitor.Combined.Cluster")]

colnames(meta) <- str_remove_all(colnames(meta), ".x")
names(meta)[names(meta) == "Age_.Carnegie_Stage." ] <- "Carnegie_Stage" 
rownames(meta) <- meta$Cell

eze <- CreateSeuratObject(counts = mat, project = "Eze_2021", meta.data=meta, assay = "RNA")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
eze[["percent.mt"]] <- PercentageFeatureSet(eze, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(eze, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(eze, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(eze, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
eze <- subset(eze, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
eze <- NormalizeData(eze)
eze <- FindVariableFeatures(eze, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(eze), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(eze)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(eze)
eze <- ScaleData(eze, features = all.genes)
eze <- RunPCA(eze, features = VariableFeatures(object = eze))
# Examine and visualize PCA results a few different ways
print(eze[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(eze, dims = 1:2, reduction = "pca")
DimPlot(eze, reduction = "pca") 
DimHeatmap(eze, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
eze <- JackStraw(eze, num.replicate = 100)
eze <- ScoreJackStraw(eze, dims = 1:20)
JackStrawPlot(eze, dims = 1:15)
ElbowPlot(eze)
eze <- FindNeighbors(eze, dims = 1:12)
eze <- FindClusters(eze, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(eze), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
eze <- RunUMAP(eze, dims = 1:12)
DimPlot(eze, reduction = "umap", group.by = "Individual")
VlnPlot(eze, features = "XIST", group.by = "Individual")
ggsave(paste0(main, "Eze_2021_XIST.pdf"))

eze@meta.data$sex <- rep("no_data", length(eze@meta.data$Cell))

eze_sex <- c("CS22"="F",    
             "CS22_2" ="M",
             "CS20"   ="F",
             "CS19"   ="F",
             "CS15_2" ="M",
             "CS14"  ="F", 
             "CS12"  ="M", 
             "CS13" ="F",  
             "CS14_3"="M")

for (id in names(eze_sex)) {
  eze@meta.data[which(eze@meta.data$Individual==id), "sex"] <- eze_sex[id]
}

VlnPlot(eze, features = "XIST", group.by = "sex")

eze@meta.data$sex_age <- paste(eze@meta.data$sex, eze@meta.data$Age, sep="_")
eze@meta.data$Individual <- str_replace_all(eze@meta.data$Individual, "_", "-")
eze@meta.data$proj <-  rep("Eze_2021", length(eze@meta.data$Cell))
eze@meta.data$id_sex_age <- paste(eze@meta.data$Individual, eze@meta.data$sex, eze@meta.data$Age, sep="_")

num_cells <- as.data.frame(table(eze$id_sex_age))
num_cells <- separate(num_cells, Var1, into=c("id", "sex", "age"), sep="_")
num_cells <- cbind("proj" = rep("Eze_2021", nrow(num_cells)), num_cells)
num_cells$age <- paste0("GW", num_cells$age)

saveRDS(eze, "UCSC/Seurat_UCSC/Eze_2021.rds")


####################################################################################################
#
# NOWAKOWSKI 2017
#
####################################################################################################


# from this tutorial UCSC CellBrowser https://cellbrowser.readthedocs.io/en/master/load.html
mat <- fread("UCSC/UCSC_downloads/exprMatrix_Nowakowski_2017.tsv.gz")
#mat <- mat[-c(which(grepl("^$", mat[,1][[1]])))]
genes <- mat[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat <- data.frame(mat[,-1], row.names=genes)

meta <- read.csv("UCSC/UCSC_downloads/new_meta_Nowakowski_2017.csv", header=T, sep=",", as.is=T, row.names=1)

rownames(meta) <- meta$Cell

nowa <- CreateSeuratObject(counts = mat, project = "Nowakowski_2017", meta.data=meta, assay = "RNA")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
nowa[["percent.mt"]] <- PercentageFeatureSet(nowa, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(nowa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(nowa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(nowa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
nowa <- subset(nowa, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
nowa <- NormalizeData(nowa)
nowa <- FindVariableFeatures(nowa, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nowa), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nowa)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(nowa)
nowa <- ScaleData(nowa, features = all.genes)
nowa <- RunPCA(nowa, features = VariableFeatures(object = nowa))
# Examine and visualize PCA results a few different ways
print(nowa[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(nowa, dims = 1:2, reduction = "pca")
DimPlot(nowa, reduction = "pca") 
DimHeatmap(nowa, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
nowa <- JackStraw(nowa, num.replicate = 100)
nowa <- ScoreJackStraw(nowa, dims = 1:20)
JackStrawPlot(nowa, dims = 1:15)
ElbowPlot(nowa)
nowa <- FindNeighbors(nowa, dims = 1:12)
nowa <- FindClusters(nowa, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(nowa), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nowa <- RunUMAP(nowa, dims = 1:12)
DimPlot(nowa, reduction = "umap")
VlnPlot(nowa, features = "XIST", group.by = "Name") + NoLegend()
ggsave(paste0(main, nowa@project.name, "_XIST.pdf"))

nowa@meta.data$sex <- rep("no_data", nrow(nowa@meta.data))

nowa_sex <- c("Sample10" = "F", 
              "Sample11"  = "M", 
              "Sample12"   = "M", 
              "Sample13"  = "F", 
              "Sample14"  = "F", 
              "Sample15"  = "F", 
              "Sample16"  = "F", 
              "Sample17"   = "M", 
              "Sample18"  = "F", 
              "Sample19"  = "F", 
              "Sample20"  = "F", 
              "Sample21" = "F", 
              "Sample22"   = "M", 
              "Sample23"  = "F", 
              "Sample24"  = "F", 
              "Sample25"   = "M", 
              "Sample26"  = "F", 
              "Sample27"   = "M", 
              "Sample28"   = "M", 
              "Sample29"  = "F", 
              "Sample3"   = "F", 
              "Sample30"   = "M", 
              "Sample31"   = "M", 
              "Sample32"  = "M", 
              "Sample33"   = "M", 
              "Sample34"   = "M", 
              "Sample35"   = "M", 
              "Sample37"   = "M", 
              "Sample38"  = "F", 
              "Sample39"   = "M", 
              "Sample4"   = "F", 
              "Sample40"   = "M", 
              "Sample41"   = "M", 
              "Sample42"  = "F", 
              "Sample43"   = "M", 
              "Sample44" = "F", 
              "Sample45"   = "M", 
              "Sample46"   = "M", 
              "Sample47"  = "M", 
              "Sample48"  = "F", 
              "Sample5"    = "M", 
              "Sample6"    = "M", 
              "Sample7"   = "F", 
              "Sample8"   = "F", 
              "Sample9"   = "M"
)

for (id in names(nowa_sex)) {
  nowa@meta.data[which(nowa@meta.data$Name==id), "sex"] <- nowa_sex[id]
}

VlnPlot(nowa, features = "XIST", group.by = "sex")

nowa@meta.data$sex_age <- paste(nowa@meta.data$sex, nowa@meta.data$Age_in_Weeks, sep="_")
nowa@meta.data$proj <-  rep(nowa@project.name, nrow(nowa@meta.data))
nowa@meta.data$id_sex_age <- paste(nowa@meta.data$Name, nowa@meta.data$sex, nowa@meta.data$Age_in_Weeks, sep="_")


num_cells3 <- as.data.frame(table(nowa$id_sex_age))
num_cells3 <- separate(num_cells3, Var1, into=c("id", "sex", "age"), sep="_")
num_cells3 <- cbind("proj" = rep(nowa@project.name, nrow(num_cells3)), num_cells3)
num_cells3$age <- paste0("GW", num_cells3$age)

num_cells <- read.csv("UCSC/outputs/num_cells_per_project.csv")
num_cells[,1] <- NULL
num_cells <- rbind(num_cells, num_cells3)
write.csv(num_cells, "UCSC/outputs/num_cells_per_project.csv")

saveRDS(nowa, paste0("UCSC/Seurat_UCSC/", nowa@project.name, ".rds"))

####################################################################################################
#
# PLOTS
#
####################################################################################################

num_cells <- read.csv("UCSC/outputs/num_cells_per_project.csv")
num_cells[,1] <- NULL
num_cells$age <- as.factor(as.numeric(str_remove_all(num_cells$age, "GW")))

p1 <- ggplot(num_cells, aes(age, fill=sex)) +
  geom_bar() +
  scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="GW Age", y="Number of samples", fill="Sex") +
  facet_wrap(~proj, scales = "free", nrow=1, drop = TRUE) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

p2 <- ggplot(num_cells, aes(age, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  #scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="GW Age", y="Number of cells/sample", fill="Sex") +
  facet_wrap(~proj, scales = "free", nrow=1) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

num_cells_proj <- ggarrange(p1, p2, common.legend = T, legend = "bottom", nrow = 2)

pdf(paste0(main, "num_samples_and_cells_per_proj.pdf"), width = 10)
print(num_cells_proj)
dev.off()