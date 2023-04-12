library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(dplyr)
`%!in%` <- Negate(`%in%`)

rds_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/"
main_deg <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/"
main_scenic <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/"

main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/"

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


# from this tutorial UCSC CellBrowser https://cellbrowser.readthedocs.io/en/master/load.html

####################################################################################################
#
# VELMESHEV 2022
#
####################################################################################################

####  Separate cell barcodes by age group - prepare indexes for file splitting

meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
split_barcodes <- list()
for (age in unique(meta$age)) {
  barcodes_id <- meta[which(meta$age==age), "cell"]
  #split_barcodes <- append(split_barcodes, list((paste(c(age, "gene", barco(des_id), collapse = ","))))
  split_barcodes <- append(split_barcodes, list((paste(c("gene", barcodes_id), collapse = ","))))
}
Reduce(intersect, split_barcodes)

sink("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_barcodes.txt")
print(split_barcodes)
sink()


names(split_barcodes) <- str_replace_all(unique(meta$age), c(" "="_", "-"="_"))


meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
split_index <- list()
index_names <- c()
problematic <- c("0-1 years", "3rd trimester", "2nd trimester")
for (age in unique(meta$age)) {
  indexes <- sort(as.numeric(rownames(meta[which(meta$age==age),])))
  if (age %in% problematic) {
    indexes <- indexes + 1
    half <- length(indexes)/2
    if (half %% 1 == 0) {
      ind1 <- indexes[1:half]
      ind2 <- indexes[(half+1):(length(indexes))]
    } else {
      ind1 <- indexes[1:(half - 0.5)]
      ind2 <- indexes[(half + 0.5):(length(indexes))]
    }
    split_index <- append(split_index, list(c(1,ind1)))
    split_index <- append(split_index, list(c(1,ind2)))
    ind1_name <- paste0(str_replace_all(age, c(" "="_", "-"="_")), "_1")
    index_names <- c(index_names, ind1_name)
    ind2_name <- paste0(str_replace_all(age, c(" "="_", "-"="_")), "_2")
    index_names <- c(index_names, ind2_name)
  } else {
    #split_index <- append(split_index, list(paste((rownames(meta[which(meta$age==age),])), collapse = ",")))
    indexes <- indexes + 1
    split_index <- append(split_index, list(c(1,indexes)))
    index_names <- c(index_names, str_replace_all(age, c(" "="_", "-"="_")))
  }
}
names(split_index) <- index_names
rm(age, half, ind1, ind1_name, ind2, ind2_name, index_names, indexes)

Reduce(intersect, split_index)

lapply(1:length(names(split_index)), function(x) write.table(
  t(split_index[[x]]),
  paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_split/", names(split_index)[x], ".txt"),
  sep=",",
  eol="",
  row.names = F,
  col.names = F))

probs_files <- list.files("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_problematic", pattern = "*.txt", full.names = T)
probs_names <- list.files("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_problematic", pattern = "*.txt", full.names = F)
probs_names <- str_remove_all(probs_names, ".txt")
mat_filt <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/exprMtx_filt_Velmeshev_2022.tsv.gz")

for (i in 1:length(probs_files)) {
  prob <- read.table(probs_files[i], sep=",", header = F)
  prob <- as.numeric(prob)
  mat_prob <- mat_filt[,..prob]
  gz_prob <- gzfile(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/", probs_names[i], ".tsv.gz"), "w")
  write.table(mat_prob, gz_prob, sep="\t")
  close(gz_prob)
}


####################################################################################################
#
# ADULT
#
####################################################################################################


meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell
#meta$samples <- gsub(".*-","",meta$cell)
#meta$samples <- str_replace_all(meta$samples, "_", "-")
#write.csv(meta, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv")


meta_ad <- subset(meta, age=="Adult")

mat_ad <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/Adult.tsv.gz")
genes <- mat_ad[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_ad <- data.frame(mat_ad[,-1], row.names=genes)
colnames(mat_ad) <- rownames(meta_ad)

velm_ad <- CreateSeuratObject(counts = mat_ad, project = "Velmeshev_2022_Adult", meta.data=meta_ad, assay = "RNA")

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_ad[["percent.mt"]] <- PercentageFeatureSet(velm_ad, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_ad, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_ad, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_ad, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_ad <- subset(velm_ad, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_ad <- NormalizeData(velm_ad)
velm_ad <- FindVariableFeatures(velm_ad, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_ad), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_ad)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_ad)
velm_ad <- ScaleData(velm_ad, features = all.genes)
velm_ad <- RunPCA(velm_ad, features = VariableFeatures(object = velm_ad))
# Examine and visualize PCA results a few different ways
print(velm_ad[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_ad, dims = 1:2, reduction = "pca")
DimPlot(velm_ad, reduction = "pca") 
DimHeatmap(velm_ad, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_ad <- JackStraw(velm_ad, num.replicate = 100)
velm_ad <- ScoreJackStraw(velm_ad, dims = 1:20)
JackStrawPlot(velm_ad, dims = 1:15)
ElbowPlot(velm_ad)
# based on rprevious plots, decide the number of dimensions
velm_ad <- FindNeighbors(velm_ad, dims = 1:12)
velm_ad <- FindClusters(velm_ad, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_ad), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_ad <- RunUMAP(velm_ad, dims = 1:12)
DimPlot(velm_ad, reduction = "umap")

velm_ad@meta.data$sex_age <- paste(velm_ad@meta.data$sex, velm_ad@meta.data$age, sep="_")
velm_ad@meta.data$proj <-  rep(velm_ad@project.name, nrow(velm_ad@meta.data))
velm_ad@meta.data$id_sex_age <- paste(velm_ad@meta.data$samples, velm_ad@meta.data$sex, velm_ad@meta.data$age, sep="_")

VlnPlot(velm_ad, features = "XIST", group.by = "id_sex_age") + NoLegend()
ggsave(paste0(main, velm_ad@project.name, "_XIST.pdf"))

saveRDS(velm_ad, paste0(rds_path, velm_ad@project.name, ".rds"))

velm_num_cells <- as.data.frame(table(velm_ad$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_ad@project.name, nrow(velm_num_cells)), velm_num_cells)

write.csv(velm_num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

################

velm_ad <- readRDS(paste0(rds_path, "Velmeshev_2022_Adult.rds"))

velm_ad@meta.data$cluster_final <- rep("no_data", nrow(velm_ad@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_ad@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_ad@meta.data[which(velm_ad@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_ad, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_ad@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_ad, reduction = "umap", group.by = "cluster_final"))
dev.off()


saveRDS(velm_ad, paste0(rds_path, velm_ad@project.name, ".rds"))

################

velm_ad <- readRDS(paste0(rds_path, "Velmeshev_2022_Adult.rds"))
velm_ad@meta.data$sex_ct <- paste(velm_ad@meta.data$sex, velm_ad@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_ad$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_ad@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

write.csv(velm_num_sex_ct, paste0(main, "Velmeshev_num_sex_ct_per_age.csv"))

saveRDS(velm_ad, paste0(rds_path, velm_ad@project.name, ".rds"))

################ FOR DEGs ANALYSIS

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/RDS_preparation/00_DEGs_and_SCENIC_files_prep.R")

velm_ad <- readRDS(paste0(rds_path, "Velmeshev_2022_Adult.rds"))

num_sex_ct <- NumSexCt(velm_ad)

min_num_cells <- c(10,50,100)

velm_ad_deg <- paste0(main_deg, velm_ad@project.name, "/outputs/")
dir.create(velm_ad_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_ad_deg)

saveRDS(velm_ad, paste0(rds_path, velm_ad@project.name, ".rds"))

############ For 02C_Conservation

velm_ad <- readRDS(paste0(rds_path, "Velmeshev_2022_Adult.rds"))

cell_info <- CreateCellInfo(velm_ad, main_deg)

expr_mat_all <- GetAssayData(velm_ad[["RNA"]], slot="data")
expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

GenerateTotGenesDf(df_list, main_deg, velm_ad)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_ad <- readRDS(paste0(rds_path, "Velmeshev_2022_Adult.rds"))

velm_ad@meta.data$sample_id <- rownames(velm_ad@meta.data)
velm_ad@meta.data$sex_ct_sample <- paste(velm_ad@meta.data$sex_ct, velm_ad@meta.data$sample_id, sep="_")
Idents(velm_ad) <- "sex_ct_sample"

velm_ad_scenic <- paste0(main_scenic, velm_ad@project.name, "/")
dir.create(velm_ad_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_ad_scenic, velm_ad)

##### Map the samples back to the groups they belong to

expr_mat_all <- readRDS(paste0(velm_ad_scenic, "top_2000_SD_expr_matrix_",  velm_ad@project.name, ".rds"))

group_list <- list()
group_list_n <- vector()
for (group_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==group_id), "cell_id"])
  df_og_group <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  group_list <- append(group_list, list(df_og_group))
  group_list_n <-  c(group_list_n, group_id)
}
names(group_list) <- group_list_n

PlotGroupNumbers(group_list, 100, velm_ad_scenic)
PlotGroupNumbers(group_list, 500, velm_ad_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_ad_scenic, "top_2000_SD_expr_matrix_",  velm_ad@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_ad@project.name, "/cell_info_", velm_ad@project.name, ".csv"))
cell_info$X <- NULL

df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_ad_scenic)


####################################################################################################
#
# 3RD TRIMESTER
#
####################################################################################################


meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell

meta_3rd_trim <- subset(meta, age=="3rd trimester")

mat_3rd_trim_1 <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/3rd_trimester_1.tsv.gz")
mat_3rd_trim_2 <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/3rd_trimester_2.tsv.gz")
mat_3rd_trim_2[,1] <- NULL
colnames(mat_3rd_trim_1) <- c("gene", (seq(1,(ncol(mat_3rd_trim_1)-1))))
colnames(mat_3rd_trim_2) <- c("gene", (seq(ncol(mat_3rd_trim_1), (ncol(mat_3rd_trim_1) + ncol(mat_3rd_trim_2) - 2))))
mat_3rd_trim <- merge(mat_3rd_trim_1, mat_3rd_trim_2, by="gene")
rm(mat_3rd_trim_1, mat_3rd_trim_2)

genes <- mat_3rd_trim[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_3rd_trim <- data.frame(mat_3rd_trim[,-1], row.names=genes)
colnames(mat_3rd_trim) <- rownames(meta_3rd_trim)

velm_3rd_trim <- CreateSeuratObject(counts = mat_3rd_trim, project = "Velmeshev_2022_3rd_trimester", meta.data=meta_3rd_trim, assay = "RNA")

saveRDS(velm_3rd_trim, paste0(rds_path, velm_3rd_trim@project.name, ".rds"))

rm(meta_3rd_trim, mat_3rd_trim)


# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_3rd_trim[["percent.mt"]] <- PercentageFeatureSet(velm_3rd_trim, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_3rd_trim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_3rd_trim, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_3rd_trim, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_3rd_trim <- subset(velm_3rd_trim, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_3rd_trim <- NormalizeData(velm_3rd_trim)
velm_3rd_trim <- FindVariableFeatures(velm_3rd_trim, selection.method = "vst", nfeatures = 2000)
saveRDS(velm_3rd_trim, paste0(rds_path, velm_3rd_trim@project.name, ".rds"))
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_3rd_trim), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_3rd_trim)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_3rd_trim)
velm_3rd_trim <- ScaleData(velm_3rd_trim, features = all.genes)
velm_3rd_trim <- RunPCA(velm_3rd_trim, features = VariableFeatures(object = velm_3rd_trim))
# Examine and visualize PCA results a few different ways
print(velm_3rd_trim[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_3rd_trim, dims = 1:2, reduction = "pca")
DimPlot(velm_3rd_trim, reduction = "pca") 
DimHeatmap(velm_3rd_trim, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_3rd_trim <- JackStraw(velm_3rd_trim, num.replicate = 100)
velm_3rd_trim <- ScoreJackStraw(velm_3rd_trim, dims = 1:20)
JackStrawPlot(velm_3rd_trim, dims = 1:15)
ElbowPlot(velm_3rd_trim)
#saveRDS(velm_3rd_trim, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_3rd_trim@project.name, ".rds"))
# based on rprevious plots, decide the number of dimensions
velm_3rd_trim <- FindNeighbors(velm_3rd_trim, dims = 1:12)
velm_3rd_trim <- FindClusters(velm_3rd_trim, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_3rd_trim), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_3rd_trim <- RunUMAP(velm_3rd_trim, dims = 1:12)
DimPlot(velm_3rd_trim, reduction = "umap")

velm_3rd_trim@meta.data$sex_age <- paste(velm_3rd_trim@meta.data$sex, velm_3rd_trim@meta.data$age, sep="_")
velm_3rd_trim@meta.data$proj <-  rep(velm_3rd_trim@project.name, nrow(velm_3rd_trim@meta.data))
velm_3rd_trim@meta.data$id_sex_age <- paste(velm_3rd_trim@meta.data$samples, velm_3rd_trim@meta.data$sex, velm_3rd_trim@meta.data$age, sep="_")

VlnPlot(velm_3rd_trim, features = "XIST", group.by = "id_sex_age")
ggsave(paste0(main, velm_3rd_trim@project.name, "_XIST.pdf"))

saveRDS(velm_3rd_trim, paste0(rds_path, velm_3rd_trim@project.name, ".rds"))

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

velm_num_cells <- as.data.frame(table(velm_3rd_trim$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_3rd_trim@project.name, nrow(velm_num_cells)), velm_num_cells)

num_cells <- rbind(num_cells, velm_num_cells)

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

################

velm_3rd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_3rd_trimester.rds"))

velm_3rd_trim@meta.data$cluster_final <- rep("no_data", nrow(velm_3rd_trim@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_3rd_trim@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_3rd_trim@meta.data[which(velm_3rd_trim@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_3rd_trim, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_3rd_trim@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_3rd_trim, reduction = "umap", group.by = "cluster_final"))
dev.off()


saveRDS(velm_3rd_trim, paste0(rds_path, velm_3rd_trim@project.name, ".rds"))

################

velm_3rd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_3rd_trimester.rds"))

velm_3rd_trim@meta.data$sex_ct <- paste(velm_3rd_trim@meta.data$sex, velm_3rd_trim@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_3rd_trim$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_3rd_trim@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

velm_num_sex_ct_all <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")
velm_num_sex_ct_all[,1] <- NULL
velm_num_sex_ct_all <- rbind(velm_num_sex_ct_all, velm_num_sex_ct)

write.csv(velm_num_sex_ct_all, paste0(main, "Velmeshev_num_sex_ct_per_age.csv"))


####################################################################################################
#
# 0-1 YEARS
#
####################################################################################################

#meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
#rownames(meta) <- meta$cell

meta_1st_year <- subset(meta, age=="0-1 years")

mat_1st_year_1 <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/0_1_years_1.tsv.gz")
mat_1st_year_2 <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/0_1_years_2.tsv.gz")
mat_1st_year_1[,1] <- NULL
colnames(mat_1st_year_1) <- c("gene", (seq(1,(ncol(mat_1st_year_1)-1))))
colnames(mat_1st_year_2) <- c("gene", (seq(ncol(mat_1st_year_1), (ncol(mat_1st_year_1) + ncol(mat_1st_year_2) - 2))))
mat_1st_year <- merge(mat_1st_year_1, mat_1st_year_2, by="gene")
rm(mat_1st_year_1, mat_1st_year_2)

genes <- mat_1st_year[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_1st_year <- data.frame(mat_1st_year[,-1], row.names=genes)
colnames(mat_1st_year) <- rownames(meta_1st_year)

velm_1st_year <- CreateSeuratObject(counts = mat_1st_year, project = "Velmeshev_2022_0_1_years", meta.data=meta_1st_year, assay = "RNA")

rm(meta_1st_year, mat_1st_year)

saveRDS(velm_1st_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_1st_year@project.name, ".rds"))


# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_1st_year[["percent.mt"]] <- PercentageFeatureSet(velm_1st_year, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_1st_year, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_1st_year, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_1st_year, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_1st_year <- subset(velm_1st_year, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_1st_year <- NormalizeData(velm_1st_year)
velm_1st_year <- FindVariableFeatures(velm_1st_year, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_1st_year), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_1st_year)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_1st_year)
velm_1st_year <- ScaleData(velm_1st_year, features = all.genes)
velm_1st_year <- RunPCA(velm_1st_year, features = VariableFeatures(object = velm_1st_year))
# Examine and visualize PCA results a few different ways
print(velm_1st_year[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_1st_year, dims = 1:2, reduction = "pca")
DimPlot(velm_1st_year, reduction = "pca") 
DimHeatmap(velm_1st_year, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_1st_year <- JackStraw(velm_1st_year, num.replicate = 100)
velm_1st_year <- ScoreJackStraw(velm_1st_year, dims = 1:20)
JackStrawPlot(velm_1st_year, dims = 1:15)
ElbowPlot(velm_1st_year)
# based on rprevious plots, decide the number of dimensions
velm_1st_year <- FindNeighbors(velm_1st_year, dims = 1:13)
velm_1st_year <- FindClusters(velm_1st_year, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_1st_year), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_1st_year <- RunUMAP(velm_1st_year, dims = 1:12)
DimPlot(velm_1st_year, reduction = "umap")

velm_1st_year@meta.data$sex_age <- paste(velm_1st_year@meta.data$sex, velm_1st_year@meta.data$age, sep="_")
velm_1st_year@meta.data$proj <-  rep(velm_1st_year@project.name, nrow(velm_1st_year@meta.data))
velm_1st_year@meta.data$id_sex_age <- paste(velm_1st_year@meta.data$samples, velm_1st_year@meta.data$sex, velm_1st_year@meta.data$age, sep="_")

VlnPlot(velm_1st_year, features = "XIST", group.by = "id_sex_age") + NoLegend()
ggsave(paste0(main, velm_1st_year@project.name, "_XIST.pdf"))

saveRDS(velm_1st_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_1st_year@project.name, ".rds"))

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

velm_num_cells <- as.data.frame(table(velm_1st_year$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_1st_year@project.name, nrow(velm_num_cells)), velm_num_cells)

num_cells <- rbind(num_cells, velm_num_cells)

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

################


velm_1st_year <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_0_1_years.rds")

velm_1st_year@meta.data$cluster_final <- rep("no_data", nrow(velm_1st_year@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_1st_year@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_1st_year@meta.data[which(velm_1st_year@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_1st_year, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_1st_year@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_1st_year, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_1st_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_1st_year@project.name, ".rds"))

################

velm_1st_year <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_0_1_years.rds")

velm_1st_year@meta.data$sex_ct <- paste(velm_1st_year@meta.data$sex, velm_1st_year@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_1st_year$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_1st_year@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

velm_num_sex_ct_all <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")
velm_num_sex_ct_all[,1] <- NULL
velm_num_sex_ct_all <- rbind(velm_num_sex_ct_all, velm_num_sex_ct)

write.csv(velm_num_sex_ct_all, paste0(main, "Velmeshev_num_sex_ct_per_age.csv"))

################ FOR DEGs ANALYSIS

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/RDS_preparation/00_DEGs_and_SCENIC_files_prep.R")

velm_1st_year <- readRDS(paste0(rds_path, "Velmeshev_2022_0_1_years.rds"))

velm_1st_year@meta.data$sex_ct <- paste(velm_1st_year@meta.data$sex, velm_1st_year@meta.data$cluster_final, sep="_")


num_sex_ct <- NumSexCt(velm_1st_year)

min_num_cells <- c(10,50,100)

velm_1st_year_deg <- paste0(main_deg, velm_1st_year@project.name, "/outputs/")
dir.create(velm_1st_year_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_1st_year_deg)

############ For 02C_Conservation

cell_info <- CreateCellInfo(velm_1st_year, main_deg)

expr_mat_all <- GetAssayData(velm_1st_year[["RNA"]], slot="data")
expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

GenerateTotGenesDf(df_list, main_deg, velm_1st_year)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_1st_year@meta.data$sample_id <- rownames(velm_1st_year@meta.data)
velm_1st_year@meta.data$sex_ct_sample <- paste(velm_1st_year@meta.data$sex_ct, velm_1st_year@meta.data$sample_id, sep="_")
Idents(velm_1st_year) <- "sex_ct_sample"

velm_1st_year_scenic <- paste0(main_scenic, velm_1st_year@project.name, "/")
dir.create(velm_1st_year_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_1st_year_scenic, velm_1st_year)

##### Map the samples back to the groups they belong to

group_list <- list()
group_list_n <- vector()
for (group_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==group_id), "cell_id"])
  df_og_group <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  group_list <- append(group_list, list(df_og_group))
  group_list_n <-  c(group_list_n, group_id)
}
names(group_list) <- group_list_n

PlotGroupNumbers(group_list, 100, velm_1st_year_scenic)
PlotGroupNumbers(group_list, 500, velm_1st_year_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_1st_year_scenic, "top_2000_SD_expr_matrix_",  velm_1st_year@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_1st_year@project.name, "/cell_info_", velm_1st_year@project.name, ".csv"))
cell_info$X <- NULL

df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_1st_year_scenic)

saveRDS(velm_1st_year, paste0(rds_path, velm_1st_year@project.name, ".rds"))


####################################################################################################
#
# 1-2- YEARS
#
####################################################################################################


#meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
#rownames(meta) <- meta$cell

meta_2nd_year <- subset(meta, age=="1-2 years")

mat_2nd_year <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/1_2_years.tsv.gz")
genes <- mat_2nd_year[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_2nd_year <- data.frame(mat_2nd_year[,-1], row.names=genes)
colnames(mat_2nd_year) <- rownames(meta_2nd_year)

velm_2nd_year <- CreateSeuratObject(counts = mat_2nd_year, project = "Velmeshev_2022_1_2_years", meta.data=meta_2nd_year, assay = "RNA")

rm(meta_2nd_year, mat_2nd_year)

saveRDS(velm_2nd_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2nd_year@project.name, ".rds"))

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_2nd_year[["percent.mt"]] <- PercentageFeatureSet(velm_2nd_year, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_2nd_year, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_2nd_year, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_2nd_year, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_2nd_year <- subset(velm_2nd_year, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_2nd_year <- NormalizeData(velm_2nd_year)
velm_2nd_year <- FindVariableFeatures(velm_2nd_year, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_2nd_year), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_2nd_year)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_2nd_year)
velm_2nd_year <- ScaleData(velm_2nd_year, features = all.genes)
velm_2nd_year <- RunPCA(velm_2nd_year, features = VariableFeatures(object = velm_2nd_year))
# Examine and visualize PCA results a few different ways
print(velm_2nd_year[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_2nd_year, dims = 1:2, reduction = "pca")
DimPlot(velm_2nd_year, reduction = "pca") 
DimHeatmap(velm_2nd_year, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_2nd_year <- JackStraw(velm_2nd_year, num.replicate = 100)
velm_2nd_year <- ScoreJackStraw(velm_2nd_year, dims = 1:20)
JackStrawPlot(velm_2nd_year, dims = 1:15)
ElbowPlot(velm_2nd_year)
# based on rprevious plots, decide the number of dimensions
velm_2nd_year <- FindNeighbors(velm_2nd_year, dims = 1:12)
velm_2nd_year <- FindClusters(velm_2nd_year, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_2nd_year), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_2nd_year <- RunUMAP(velm_2nd_year, dims = 1:12)
DimPlot(velm_2nd_year, reduction = "umap")

velm_2nd_year@meta.data$sex_age <- paste(velm_2nd_year@meta.data$sex, velm_2nd_year@meta.data$age, sep="_")
velm_2nd_year@meta.data$proj <-  rep(velm_2nd_year@project.name, nrow(velm_2nd_year@meta.data))
velm_2nd_year@meta.data$id_sex_age <- paste(velm_2nd_year@meta.data$samples, velm_2nd_year@meta.data$sex, velm_2nd_year@meta.data$age, sep="_")

VlnPlot(velm_2nd_year, features = "XIST", group.by = "id_sex_age") + NoLegend()
ggsave(paste0(main, velm_2nd_year@project.name, "_XIST.pdf"))

saveRDS(velm_2nd_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2nd_year@project.name, ".rds"))

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

velm_num_cells <- as.data.frame(table(velm_2nd_year$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_2nd_year@project.name, nrow(velm_num_cells)), velm_num_cells)

num_cells <- rbind(num_cells, velm_num_cells)

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

################

velm_2nd_year <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_1_2_years.rds")

velm_2nd_year@meta.data$cluster_final <- rep("no_data", nrow(velm_2nd_year@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_2nd_year@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_2nd_year@meta.data[which(velm_2nd_year@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_2nd_year, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_2nd_year@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_2nd_year, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_2nd_year, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2nd_year@project.name, ".rds"))

################

velm_2nd_year <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_1_2_years.rds")

velm_2nd_year@meta.data$sex_ct <- paste(velm_2nd_year@meta.data$sex, velm_2nd_year@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_2nd_year$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_2nd_year@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

velm_num_sex_ct_all <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")
velm_num_sex_ct_all[,1] <- NULL
velm_num_sex_ct_all <- rbind(velm_num_sex_ct_all, velm_num_sex_ct)

write.csv(velm_num_sex_ct_all, paste0(main, "Velmeshev_num_sex_ct_per_age.csv"))

################ FOR DEGs ANALYSIS

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/RDS_preparation/00_DEGs_and_SCENIC_files_prep.R")

velm_2nd_year <- readRDS(paste0(rds_path, "Velmeshev_2022_1_2_years.rds"))

velm_2nd_year@meta.data$sex_ct <- paste(velm_2nd_year@meta.data$sex, velm_2nd_year@meta.data$cluster_final, sep="_")

num_sex_ct <- NumSexCt(velm_2nd_year)

min_num_cells <- c(10,50,100)

velm_2nd_year_deg <- paste0(main_deg, velm_2nd_year@project.name, "/outputs/")
dir.create(velm_2nd_year_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_2nd_year_deg)


############ For 02C_Conservation

cell_info <- CreateCellInfo(velm_2nd_year, main_deg)

expr_mat_all <- GetAssayData(velm_2nd_year[["RNA"]], slot="data")
expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

GenerateTotGenesDf(df_list, main_deg, velm_2nd_year)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_2nd_year@meta.data$sample_id <- rownames(velm_2nd_year@meta.data)
velm_2nd_year@meta.data$sex_ct_sample <- paste(velm_2nd_year@meta.data$sex_ct, velm_2nd_year@meta.data$sample_id, sep="_")
Idents(velm_2nd_year) <- "sex_ct_sample"

velm_2nd_year_scenic <- paste0(main_scenic, velm_2nd_year@project.name, "/")
dir.create(velm_2nd_year_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_2nd_year_scenic, velm_2nd_year)

##### Map the samples back to the groups they belong to

group_list <- list()
group_list_n <- vector()
for (group_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==group_id), "cell_id"])
  df_og_group <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  group_list <- append(group_list, list(df_og_group))
  group_list_n <-  c(group_list_n, group_id)
}
names(group_list) <- group_list_n

PlotGroupNumbers(group_list, 100, velm_2nd_year_scenic)
PlotGroupNumbers(group_list, 500, velm_2nd_year_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_2nd_year_scenic, "top_2000_SD_expr_matrix_",  velm_2nd_year@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_2nd_year@project.name, "/cell_info_", velm_2nd_year@project.name, ".csv"))
cell_info$X <- NULL

df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_2nd_year_scenic)

saveRDS(velm_2nd_year, paste0(rds_path, velm_2nd_year@project.name, ".rds"))

####################################################################################################
#
# 2-4 YEARS
#
####################################################################################################


meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell
meta_2_4_years <- subset(meta, age=="2-4 years")

mat_2_4_years <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/2_4_years.tsv.gz")
genes <- mat_2_4_years[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_2_4_years <- data.frame(mat_2_4_years[,-1], row.names=genes)
colnames(mat_2_4_years) <- rownames(meta_2_4_years)

velm_2_4_years <- CreateSeuratObject(counts = mat_2_4_years, project = "Velmeshev_2022_2_4_years", meta.data=meta_2_4_years, assay = "RNA")

rm(meta_2_4_years, mat_2_4_years, genes)

#saveRDS(velm_2_4_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2_4_years@project.name, ".rds"))

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_2_4_years[["percent.mt"]] <- PercentageFeatureSet(velm_2_4_years, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_2_4_years, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_2_4_years, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_2_4_years, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_2_4_years <- subset(velm_2_4_years, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_2_4_years <- NormalizeData(velm_2_4_years)
velm_2_4_years <- FindVariableFeatures(velm_2_4_years, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_2_4_years), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_2_4_years)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_2_4_years)
velm_2_4_years <- ScaleData(velm_2_4_years, features = all.genes)
velm_2_4_years <- RunPCA(velm_2_4_years, features = VariableFeatures(object = velm_2_4_years))
# Examine and visualize PCA results a few different ways
print(velm_2_4_years[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_2_4_years, dims = 1:2, reduction = "pca")
DimPlot(velm_2_4_years, reduction = "pca") 
DimHeatmap(velm_2_4_years, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_2_4_years <- JackStraw(velm_2_4_years, num.replicate = 100)
velm_2_4_years <- ScoreJackStraw(velm_2_4_years, dims = 1:20)
JackStrawPlot(velm_2_4_years, dims = 1:15)
ElbowPlot(velm_2_4_years)
# based on rprevious plots, decide the number of dimensions
velm_2_4_years <- FindNeighbors(velm_2_4_years, dims = 1:13)
velm_2_4_years <- FindClusters(velm_2_4_years, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_2_4_years), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_2_4_years <- RunUMAP(velm_2_4_years, dims = 1:12)
DimPlot(velm_2_4_years, reduction = "umap")

velm_2_4_years@meta.data$sex_age <- paste(velm_2_4_years@meta.data$sex, velm_2_4_years@meta.data$age, sep="_")
velm_2_4_years@meta.data$proj <-  rep(velm_2_4_years@project.name, nrow(velm_2_4_years@meta.data))
velm_2_4_years@meta.data$id_sex_age <- paste(velm_2_4_years@meta.data$samples, velm_2_4_years@meta.data$sex, velm_2_4_years@meta.data$age, sep="_")

VlnPlot(velm_2_4_years, features = "XIST", group.by = "id_sex_age") + NoLegend()
ggsave(paste0(main, velm_2_4_years@project.name, "_XIST.pdf"))

saveRDS(velm_2_4_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2_4_years@project.name, ".rds"))

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

velm_num_cells <- as.data.frame(table(velm_2_4_years$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_2_4_years@project.name, nrow(velm_num_cells)), velm_num_cells)

num_cells <- rbind(num_cells, velm_num_cells)

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

################

velm_2_4_years <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_2_4_years.rds")

velm_2_4_years@meta.data$cluster_final <- rep("no_data", nrow(velm_2_4_years@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_2_4_years@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_2_4_years@meta.data[which(velm_2_4_years@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_2_4_years, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_2_4_years@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_2_4_years, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_2_4_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_2_4_years@project.name, ".rds"))

################

velm_2_4_years <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Velmeshev_2022_2_4_years.rds")

velm_2_4_years@meta.data$sex_ct <- paste(velm_2_4_years@meta.data$sex, velm_2_4_years@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_2_4_years$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_2_4_years@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

velm_num_sex_ct_all <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")
velm_num_sex_ct_all[,1] <- NULL
velm_num_sex_ct_all <- rbind(velm_num_sex_ct_all, velm_num_sex_ct)

write.csv(velm_num_sex_ct_all, paste0(main, "Velmeshev_num_sex_ct_per_age.csv"))

################ FOR DEGs ANALYSIS

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/RDS_preparation/00_DEGs_and_SCENIC_files_prep.R")

velm_2_4_years <- readRDS(paste0(rds_path, "Velmeshev_2022_2_4_years.rds"))

velm_2_4_years@meta.data$sex_ct <- paste(velm_2_4_years@meta.data$sex, velm_2_4_years@meta.data$cluster_final, sep="_")

num_sex_ct <- NumSexCt(velm_2_4_years)

min_num_cells <- c(10,50,100)

velm_2_4_years_deg <- paste0(main_deg, velm_2_4_years@project.name, "/outputs/")
dir.create(velm_2_4_years_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_2_4_years_deg)

############ For 02C_Conservation

cell_info <- CreateCellInfo(velm_2_4_years, main_deg)

expr_mat_all <- GetAssayData(velm_2_4_years[["RNA"]], slot="data")
expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

GenerateTotGenesDf(df_list, main_deg, velm_2_4_years)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_2_4_years@meta.data$sample_id <- rownames(velm_2_4_years@meta.data)
velm_2_4_years@meta.data$sex_ct_sample <- paste(velm_2_4_years@meta.data$sex_ct, velm_2_4_years@meta.data$sample_id, sep="_")
Idents(velm_2_4_years) <- "sex_ct_sample"

velm_2_4_years_scenic <- paste0(main_scenic, velm_2_4_years@project.name, "/")
dir.create(velm_2_4_years_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_2_4_years_scenic, velm_2_4_years)

##### Map the samples back to the groups they belong to

group_list <- list()
group_list_n <- vector()
for (group_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==group_id), "cell_id"])
  df_og_group <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  group_list <- append(group_list, list(df_og_group))
  group_list_n <-  c(group_list_n, group_id)
}
names(group_list) <- group_list_n

PlotGroupNumbers(group_list, 100, velm_2_4_years_scenic)
PlotGroupNumbers(group_list, 500, velm_2_4_years_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_2_4_years_scenic, "top_2000_SD_expr_matrix_",  velm_2_4_years@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_2_4_years@project.name, "/cell_info_", velm_2_4_years@project.name, ".csv"))
cell_info$X <- NULL


df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_2_4_years_scenic)


saveRDS(velm_2_4_years, paste0(rds_path, velm_2_4_years@project.name, ".rds"))

####################################################################################################
#
# 4-10 YEARS
#
####################################################################################################


meta <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
rownames(meta) <- meta$cell

meta_4_10_years <- subset(meta, age=="4-10 years")

mat_4_10_years <- fread("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/Velmeshev_outs/4_10_years.tsv.gz")
genes <- mat_4_10_years[,1][[1]]
genes <- gsub(".+[|]", "", genes)
mat_4_10_years <- data.frame(mat_4_10_years[,-1], row.names=genes)
colnames(mat_4_10_years) <- rownames(meta_4_10_years)

velm_4_10_years <- CreateSeuratObject(counts = mat_4_10_years, project = "Velmeshev_2022_4_10_years", meta.data=meta_4_10_years, assay = "RNA")

rm(meta_4_10_years, mat_4_10_years, genes)

#saveRDS(velm_4_10_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_4_10_years@project.name, ".rds"))

# Seurat tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
velm_4_10_years[["percent.mt"]] <- PercentageFeatureSet(velm_4_10_years, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(velm_4_10_years, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(velm_4_10_years, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(velm_4_10_years, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
velm_4_10_years <- subset(velm_4_10_years, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
velm_4_10_years <- NormalizeData(velm_4_10_years)
velm_4_10_years <- FindVariableFeatures(velm_4_10_years, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(velm_4_10_years), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(velm_4_10_years)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(velm_4_10_years)
velm_4_10_years <- ScaleData(velm_4_10_years, features = all.genes)
velm_4_10_years <- RunPCA(velm_4_10_years, features = VariableFeatures(object = velm_4_10_years))
# Examine and visualize PCA results a few different ways
print(velm_4_10_years[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(velm_4_10_years, dims = 1:2, reduction = "pca")
DimPlot(velm_4_10_years, reduction = "pca") 
DimHeatmap(velm_4_10_years, dims = 1, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
velm_4_10_years <- JackStraw(velm_4_10_years, num.replicate = 100)
velm_4_10_years <- ScoreJackStraw(velm_4_10_years, dims = 1:20)
JackStrawPlot(velm_4_10_years, dims = 1:15)
ElbowPlot(velm_4_10_years)
# based on rprevious plots, decide the number of dimensions
velm_4_10_years <- FindNeighbors(velm_4_10_years, dims = 1:12)
velm_4_10_years <- FindClusters(velm_4_10_years, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
#head(Idents(velm_4_10_years), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
velm_4_10_years <- RunUMAP(velm_4_10_years, dims = 1:12)
DimPlot(velm_4_10_years, reduction = "umap")

velm_4_10_years@meta.data$sex_age <- paste(velm_4_10_years@meta.data$sex, velm_4_10_years@meta.data$age, sep="_")
velm_4_10_years@meta.data$proj <-  rep(velm_4_10_years@project.name, nrow(velm_4_10_years@meta.data))
velm_4_10_years@meta.data$id_sex_age <- paste(velm_4_10_years@meta.data$samples, velm_4_10_years@meta.data$sex, velm_4_10_years@meta.data$age, sep="_")

VlnPlot(velm_4_10_years, features = "XIST", group.by = "id_sex_age") + NoLegend()
ggsave(paste0(main, velm_4_10_years@project.name, "_XIST.pdf"))

saveRDS(velm_4_10_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/", velm_4_10_years@project.name, ".rds"))

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

velm_num_cells <- as.data.frame(table(velm_4_10_years$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_4_10_years@project.name, nrow(velm_num_cells)), velm_num_cells)

num_cells <- rbind(num_cells, velm_num_cells)

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")

#### no more analysis on this dataset because it has only 1F! 

velm_4_10_years <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/excluded/Velmeshev_2022_4_10_years.rds")

velm_4_10_years@meta.data$cluster_final <- rep("no_data", nrow(velm_4_10_years@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_4_10_years@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_4_10_years@meta.data[which(velm_4_10_years@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}
DimPlot(velm_4_10_years, reduction = "umap", group.by = "cluster_final")

pdf(paste0(main, velm_4_10_years@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_4_10_years, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_4_10_years, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/excluded/", velm_4_10_years@project.name, ".rds"))
