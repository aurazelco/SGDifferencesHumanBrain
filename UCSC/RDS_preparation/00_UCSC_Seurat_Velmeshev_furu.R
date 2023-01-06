library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(dplyr)
`%!in%` <- Negate(`%in%`)


out_path <- "/Home/ii/auraz/data/UCSC"
rds_path <- paste0(out_path, "/Seurat_UCSC/Velmeshev/")
main_deg <- "/Home/ii/auraz/data/UCSC/outputs/DEGs/"
main_scenic <- "/Home/ii/auraz/data/UCSC/outputs/SCENIC/"

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

output_2nd_trim <- paste0(out_path, "/outputs/Velmeshev_2nd_trimester/")
dir.create(output_2nd_trim, recursive=T, showWarnings = F)
dir.create(rds_path, recursive=T, showWarnings = F)

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

saveRDS(velm_2nd_trim, paste0(rds_path, velm_2nd_trim@project.name, ".rds"))

velm_num_cells <- as.data.frame(table(velm_2nd_trim$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_2nd_trim@project.name, nrow(velm_num_cells)), velm_num_cells)
write.csv(velm_num_cells, paste0(output_2nd_trim, velm_2nd_trim@project.name, "_num_cells.csv"))

################

velm_2nd_trim@meta.data$cluster_final <- rep("no_data", nrow(velm_2nd_trim@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_2nd_trim@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_2nd_trim@meta.data[which(velm_2nd_trim@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}

pdf(paste0(output_2nd_trim, velm_2nd_trim@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_2nd_trim, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_2nd_trim, paste0(rds_path, velm_2nd_trim@project.name, ".rds"))

################ FOR DEGs ANALYSIS

velm_2nd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_2nd_trimester.rds"))

velm_2nd_trim@meta.data$sex_ct <- paste(velm_2nd_trim@meta.data$sex, velm_2nd_trim@meta.data$cluster_final, sep="_")

num_sex_ct <- as.data.frame(table(velm_2nd_trim$sex_ct))
num_sex_ct <- separate(num_sex_ct, Var1, into = c("sex" , "ct"), sep = "_", remove = F)
names(num_sex_ct)[names(num_sex_ct) == 'Var1'] <- "idents"
names(num_sex_ct)[names(num_sex_ct) == 'Freq'] <- "count"
col_factors <- c("idents", "sex","ct")
num_sex_ct[col_factors] <- lapply(num_sex_ct[col_factors], as.factor)  

FiltDF <- function(df, min_num_cells) {
  `%!in%` <- Negate(`%in%`)
  df <- droplevels(df)
  incomplete_ct <- vector()
  for (type in levels(df$ct)) {
    if ((nrow(subset(df, subset = ct==type))%%2!=0) | (any(subset(df, subset = ct==type)[,"count"] < min_num_cells))) {
      incomplete_ct <- c(incomplete_ct, type)
    }
  }
  df_filt <- df[df$ct %!in% incomplete_ct,]
  return(df_filt)
}

min_num_cells <- c(10,50,100)

main_deg <- "/Home/ii/auraz/data/UCSC/outputs/DEGs/"

velm_2nd_trim_output <- paste0(main_deg, velm_2nd_trim@project.name, "/outputs/")

dir.create(velm_2nd_trim_output, recursive = T, showWarnings = F)

for (min_cells in min_num_cells) {
  num_filt <- FiltDF(num_sex_ct, min_cells)
  write.csv(num_filt, file = paste0(velm_2nd_trim_output, "final_filt_", min_cells, ".csv"),
            row.names = F)
  pdf(paste0(velm_2nd_trim_output, "filt_counts_", min_cells, ".pdf"), 10, 15)
  print(ggplot(num_sex_ct, aes(ct, count, fill=sex)) +
          geom_bar(stat="identity", position = "dodge") + 
          labs(x="", y="Nuclei count", fill="Sex") +
          geom_hline(yintercept = min_cells, linetype="dashed") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.title.x = element_text(size=12, face="bold", colour = "black"),
                axis.text.x = element_text(size=8, colour = "black",angle = 45, vjust = 0.5, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom"))
  dev.off()
}

saveRDS(velm_2nd_trim, paste0(rds_path, velm_2nd_trim@project.name, ".rds"))


############ For 02C_Conservation

velm_2nd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_2nd_trimester.rds"))

Idents(velm_2nd_trim) <- "sex_ct"

expr_mat_all <- GetAssayData(velm_2nd_trim[["RNA"]], slot="data")

cell_info <- data.frame()
for (i in unique(velm_2nd_trim@meta.data$sex_ct)) {
  print(i)
  cell_id <- WhichCells(velm_2nd_trim, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}

write.csv(cell_info, paste0(main_deg, velm_2nd_trim@project.name, "/cell_info_", velm_2nd_trim@project.name,   ".csv"))

rm(velm_2nd_trim)
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

FiltDf <- function(df_list_dis) {
  filt_names <- vector()
  df_dis <- list()
  for (k in names(df_list_dis)) {
    if (!is.null(ncol(df_list_dis[[k]]))) {
      df_dis <- append(df_dis, list(rownames(df_list_dis[[k]][which(rowSums(as.matrix(df_list_dis[[k]]))!=0),])))
      filt_names <- c(filt_names, k)
    }
  }
  names(df_dis) <- filt_names
  cts <- vector()
  genes <- vector()
  for (id in names(df_dis)) {
    cts <- c(cts, rep(id, length(df_dis[[id]])))
    genes <- c(genes, df_dis[[id]])
  }
  tot_genes <- data.frame(cts, genes)
  return(tot_genes)
}

tot_df <- FiltDf(df_list)

tot_df <- separate(tot_df, cts, into=c("sex", "ct"), sep ="_", remove = FALSE)
colnames(tot_df)
names(tot_df)[names(tot_df) == "cts"] <- "og"

col_factors <- c("og", "sex","ct")
tot_df[col_factors] <- lapply(tot_df[col_factors], as.factor) 

sexes <- vector()
cts <- vector()
genes <- vector()
for (sex_id in levels(tot_df$sex)) {
  for (ct_id in levels(tot_df$ct)) {
    common_genes <- tot_df[which(tot_df$sex==sex_id & tot_df$ct==ct_id), "genes"]
    genes <- c(genes, common_genes)
    sexes <- c(sexes, rep(sex_id, length(common_genes)))
    cts <- c(cts, rep(ct_id, length(common_genes)))
  }
}

tot_genes <- as.data.frame(cbind(sexes, cts, genes))
colnames(tot_genes) <- c("sex", "ct", "genes")
col_factors <- c("sex", "ct")
tot_genes[col_factors] <- lapply(tot_genes[col_factors], as.factor) 

write.csv(tot_genes, paste0(main_deg, velm_2nd_trim@project.name, "/tot_genes_ct_", velm_2nd_trim@project.name,   ".csv"))

############ For SCENIC - HIGHEST VARIABLE GENES

velm_2nd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_2nd_trimester.rds"))

velm_2nd_trim@meta.data$sample_id <- rownames(velm_2nd_trim@meta.data)

velm_2nd_trim@meta.data$sex_ct_sample <- paste(velm_2nd_trim@meta.data$sex_ct, velm_2nd_trim@meta.data$sample_id, sep="_")

Idents(velm_2nd_trim) <- "sex_ct_sample"

#expr_mat_all <- GetAssayData(velm_2nd_trim[["RNA"]], slot="data")

#cell_info <- read.csv(paste0(main_deg, velm_2nd_trim@project.name, "/cell_info_", velm_2nd_trim@project.name, ".csv"))
#cell_info$X <- NULL

#saveRDS(velm_2nd_trim, paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

#rm(velm_2nd_trim)

#expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))
expr_mat_all$SD <- rowSds(as.matrix(expr_mat_all))
expr_mat_all <- expr_mat_all[which(expr_mat_all$SD > 0), ]

if (nrow(expr_mat_all) * 0.25 > 2000) {
  expr_mat_all <- expr_mat_all[which(expr_mat_all$SD > quantile(expr_mat_all$SD)[4]), ]
} else {
  print(" less than 2k genes above third quantile")
}

# order df in descending order
expr_mat_all <- expr_mat_all[order(-expr_mat_all$SD),] 
expr_mat_all <- expr_mat_all[1:2000, ]
expr_mat_all <- cbind("Genes" = rownames(expr_mat_all), expr_mat_all)
rownames(expr_mat_all) <- NULL

expr_sums <- colSums(expr_mat_all[2:ncol(expr_mat_all)])
if (identical(length(which(expr_sums>0)), length(expr_sums))) {
  print("all columns express at least one gene")
} else {
  expr_mat_all <- expr_mat_all[ , !(names(expr_mat_all) %in% which(expr_sums>0))]
  print("calculate how many cells have been filtered out")
}

expr_mat_all <- expr_mat_all %>% 
  relocate(SD, .after = Genes)

velm_2nd_trim_scenic <- paste0("/Home/ii/auraz/data/UCSC/outputs/SCENIC/", velm_2nd_trim@project.name, "/")

dir.create(velm_2nd_trim_scenic, showWarnings = F, recursive = T)

saveRDS(expr_mat_all, paste0(velm_2nd_trim_scenic, "top_2000_SD_expr_matrix",  velm_2nd_trim@project.name, ".rds"))

##### Map the samples back to the groups they belong to

expr_mat_all <- readRDS(paste0(velm_2nd_trim_scenic, "top_2000_SD_expr_matrix_",  velm_2nd_trim@project.name, ".rds"))

group_list <- list()
group_list_n <- vector()
for (group_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==group_id), "cell_id"])
  df_og_group <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  group_list <- append(group_list, list(df_og_group))
  group_list_n <-  c(group_list_n, group_id)
}
names(group_list) <- group_list_n

remove_dfs <- function(df_list, threshold) {
  incomplete_dfs <- vector()
  for (group_id in names(df_list)) {
    if (ncol(df_list[[group_id]]) < (threshold + 1)) {
      incomplete_dfs <- c(incomplete_dfs, group_id)
    }
  }
  add_counterpart <- vector()
  for (i in incomplete_dfs) {
    if (grepl("Female", i)) {
      m_id <- str_replace(i, "Female", "Male")
      if (m_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, m_id)
      }
    } else {
      f_id <- str_replace(i, "Male", "Female")
      if (f_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, f_id)
      }
    }
  }
  incomplete_dfs <- c(incomplete_dfs, add_counterpart)
  for (i in incomplete_dfs) {
    df_list[[i]] <- NULL
  }
  return(df_list)
}

group_list100 <- remove_dfs(group_list, 100)
group_list500 <- remove_dfs(group_list, 500)

plot_group_numbers <- function(df_list, thresh) {
  ids <- as.data.frame(names(df_list))
  colnames(ids) <- c("groups")
  ids <- separate(ids, groups, into = c("sex","ct"), sep="_", remove=FALSE)
  col_factors <- c("sex","ct")
  ids[col_factors] <- lapply(ids[col_factors], as.factor)  
  ids$length_groups <- sapply(1:length(names(df_list)), function(i) ncol(df_list[[i]]))
  pdf(paste0(velm_2nd_trim_scenic, "num_filt_", thresh, "_cells.pdf"))
  print(
    ggplot(ids, aes(ct, length_groups, fill=sex)) +
      geom_bar(stat="identity", position = "dodge") + 
      geom_hline(yintercept = thresh, linetype="dashed", color = "black") +
      labs(title = paste0("Filter: ", thresh, " cells"), x="cell types", y="# of cells", fill="sex") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.5, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom",
            plot.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}
plot_group_numbers(group_list100, 100)
plot_group_numbers(group_list500, 500)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_2nd_trim_scenic, "top_2000_SD_expr_matrix_",  velm_2nd_trim@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_2nd_trim@project.name, "/cell_info_", velm_2nd_trim@project.name, ".csv"))
cell_info$X <- NULL

remove_dfs <- function(df_list, threshold) {
  incomplete_dfs <- vector()
  for (group_id in names(df_list)) {
    if (ncol(df_list[[group_id]]) < (threshold + 1)) {
      incomplete_dfs <- c(incomplete_dfs, group_id)
    }
  }
  add_counterpart <- vector()
  for (i in incomplete_dfs) {
    if (grepl("F", i)) {
      m_id <- str_replace(i, "Female", "Male")
      if (m_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, m_id)
      }
    } else {
      f_id <- str_replace(i, "Male", "Female")
      if (f_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, f_id)
      }
    }
  }
  incomplete_dfs <- c(incomplete_dfs, add_counterpart)
  for (i in incomplete_dfs) {
    df_list[[i]] <- NULL
  }
  return(df_list)
}

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all[ , (names(expr_mat_all) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

df_list100 <- df_list
df_list100 <- remove_dfs(df_list100, 100)

check_dfs <- function(group_list) {
  remove_groups <- vector()
  if (length(group_list) %% 2 != 0 ) {
    f_list <- vector()
    m_list <- vector()
    for (i in names(group_list)) {
      if (grepl("Female", i)) {
        gen <- str_remove(i, "Female_")
        f_list <- c(f_list, gen)
      } else {
        gen <- str_remove(i, "Male_")
        m_list <- c(m_list, gen)
      }
    }
    if (identical(m_list, f_list) == FALSE) {
      if (length(m_list) > length(f_list)) {
        remove_groups <- c(m_list[which(m_list %!in% f_list)], "Male")
      } else {
        remove_groups <- c(f_list[which(f_list %!in% m_list)], "Female")
      }
    }
  }
  return(remove_groups)
}

check_dfs(df_list100)

rand_sample <- function(group_list, num_sampling, num_cells, main) {
  sampled_dfs <-list()
  sampled_names <- vector()
  for (id in names(group_list)) {
    for (k in 1:num_sampling) {
      sampled <- data.frame()
      sampled <- sample(group_list[[id]][-1], num_cells)
      sampled <- cbind("Genes" = group_list[[id]]$Genes, sampled)
      sampled_dfs <- append(sampled_dfs, list(sampled))
      sampled_names <- c(sampled_names, paste(id, k, sep="_"))
    }
  }
  names(sampled_dfs) <- lapply(1:length(sampled_names), function(i) str_replace_all(sampled_names[i]," ", "_"))
  dir.create(paste0(main_scenic, "0_input_dfs/sampled_", num_cells, "_cells"), showWarnings = FALSE, recursive = T)
  lapply(1:length(names(sampled_dfs)), function(i) write.csv(sampled_dfs[[i]], 
                                                             file = paste0(main_scenic, "0_input_dfs/sampled_", num_cells, "_cells/", names(sampled_dfs)[i], ".csv"),
                                                             row.names = FALSE))
  return(sampled_dfs)
}


df_100_sampled <- rand_sample(df_list100, 3, 100, main_scenic)

####################################################################################################
#
# 10-20 YEARS
#
####################################################################################################

output_10_20_yo <- paste0(out_path, "/outputs/Velmeshev_10_20_years")
rds_path <- paste0(out_path, "/Seurat_UCSC/Velmeshev")
dir.create(output_10_20_yo, recursive=T, showWarnings = F)
dir.create(rds_path, recursive=T, showWarnings = F)

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

pdf(paste0(output_10_20_yo, velm_10_20_years@project.name,  "_XIST.pdf"))
print(VlnPlot(velm_10_20_years, features = "XIST", group.by = "id_sex_age") + NoLegend())
dev.off()

#saveRDS(velm_10_20_years, paste0(out_path, "/Seurat_UCSC/", velm_10_20_years@project.name, ".rds"))

velm_num_cells <- as.data.frame(table(velm_10_20_years$id_sex_age))
velm_num_cells <- separate(velm_num_cells, Var1, into=c("id", "sex", "age"), sep="_")
velm_num_cells <- cbind("proj" = rep(velm_10_20_years@project.name, nrow(velm_num_cells)), velm_num_cells)
write.csv(velm_num_cells, paste0(output_10_20_yo, velm_10_20_years@project.name, "_num_cells.csv"))


################

velm_10_20_years@meta.data$cluster_final <- rep("no_data", nrow(velm_10_20_years@meta.data))
present_clusters <- ann_df[which(ann_df$og_clusters %in% velm_10_20_years@meta.data$cluster),]
for (ct in unique(present_clusters$cts)) {
  print(ct)
  for (og_cl in present_clusters[which(present_clusters$cts==ct), "og_clusters"]) {
    velm_10_20_years@meta.data[which(velm_10_20_years@meta.data$cluster==og_cl), "cluster_final"] <- ct
  }
}

pdf(paste0(output_10_20_yo, velm_10_20_years@project.name,  "_cluster_final.pdf"))
print(DimPlot(velm_10_20_years, reduction = "umap", group.by = "cluster_final"))
dev.off()

saveRDS(velm_10_20_years, paste0(rds_path, velm_10_20_years@project.name, ".rds"))

################

velm_10_20_years <- readRDS(paste0(rds_path, "Velmeshev_2022_10_20_years.rds"))

velm_10_20_years@meta.data$sex_ct <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$cluster_final, sep="_")

velm_num_sex_ct <- as.data.frame(table(velm_10_20_years$sex_ct))
velm_num_sex_ct <- separate(velm_num_sex_ct, Var1, into=c("sex", "ct"), sep="_")
velm_num_sex_ct <- cbind("proj" = rep(velm_10_20_years@project.name, nrow(velm_num_sex_ct)), velm_num_sex_ct)

write.csv(velm_num_sex_ct, paste0(output_10_20_yo, velm_10_20_years@project.name, "_num_sex_ct_per_age.csv"))

################ FOR DEGs ANALYSIS

source("/Home/ii/auraz/scripts/00_DEGs_and_SCENIC_files_prep.R")

velm_10_20_years <- readRDS(paste0(rds_path, "Velmeshev_2022_10_20_years.rds"))

velm_10_20_years@meta.data$sex_ct <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$cluster_final, sep="_")

num_sex_ct <- NumSexCt(velm_10_20_years)

min_num_cells <- c(10,50,100)

velm_10_20_years_deg <- paste0(main_deg, velm_10_20_years@project.name, "/outputs/")
dir.create(velm_10_20_years_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_10_20_years_deg)

############ For 02C_Conservation

cell_info <- CreateCellInfo(velm_10_20_years, main_deg)

expr_mat_all <- GetAssayData(velm_10_20_years[["RNA"]], slot="data")
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

GenerateTotGenesDf(df_list, main_deg, velm_10_20_years)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_10_20_years@meta.data$sample_id <- rownames(velm_10_20_years@meta.data)
velm_10_20_years@meta.data$sex_ct_sample <- paste(velm_10_20_years@meta.data$sex_ct, velm_10_20_years@meta.data$sample_id, sep="_")
Idents(velm_10_20_years) <- "sex_ct_sample"

velm_10_20_years_scenic <- paste0(main_scenic, velm_10_20_years@project.name, "/")
dir.create(velm_10_20_years_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_10_20_years_scenic, velm_10_20_years)

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

PlotGroupNumbers(group_list, 100, velm_10_20_years_scenic)
PlotGroupNumbers(group_list, 500, velm_10_20_years_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_10_20_years_scenic, "top_2000_SD_expr_matrix_",  velm_10_20_years@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_10_20_years@project.name, "/cell_info_", velm_10_20_years@project.name, ".csv"))
cell_info$X <- NULL


df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_10_20_years_scenic)





####################################################################################################
#
# 3RD TRIMESTER
#
####################################################################################################

################ FOR DEGs ANALYSIS

source("/Home/ii/auraz/scripts/00_DEGs_and_SCENIC_files_prep.R")

velm_3rd_trim <- readRDS(paste0(rds_path, "Velmeshev_2022_3rd_trimester.rds"))

velm_3rd_trim@meta.data$sex_ct <- paste(velm_3rd_trim@meta.data$sex, velm_3rd_trim@meta.data$cluster_final, sep="_")

num_sex_ct <- NumSexCt(velm_3rd_trim)

min_num_cells <- c(10,50,100)

velm_3rd_trim_deg <- paste0(main_deg, velm_3rd_trim@project.name, "/outputs/")
dir.create(velm_3rd_trim_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_3rd_trim_deg)

############ For 02C_Conservation

cell_info <- CreateCellInfo(velm_3rd_trim, main_deg)

expr_mat_all <- GetAssayData(velm_3rd_trim[["RNA"]], slot="data")
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

GenerateTotGenesDf(df_list, main_deg, velm_3rd_trim)

############ For SCENIC - HIGHEST VARIABLE GENES

velm_3rd_trim@meta.data$sample_id <- rownames(velm_3rd_trim@meta.data)
velm_3rd_trim@meta.data$sex_ct_sample <- paste(velm_3rd_trim@meta.data$sex_ct, velm_3rd_trim@meta.data$sample_id, sep="_")
Idents(velm_3rd_trim) <- "sex_ct_sample"

velm_3rd_trim_scenic <- paste0(main_scenic, velm_3rd_trim@project.name, "/")
dir.create(velm_3rd_trim_scenic, recursive = T, showWarnings = F)

Top2000ExprMtx(expr_mat_all, velm_3rd_trim_scenic, velm_3rd_trim)

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

PlotGroupNumbers(group_list, 100, velm_3rd_trim_scenic)
PlotGroupNumbers(group_list, 500, velm_3rd_trim_scenic)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(velm_3rd_trim_scenic, "top_2000_SD_expr_matrix_",  velm_3rd_trim@project.name, ".rds"))
cell_info <- read.csv(paste0(main_deg, velm_3rd_trim@project.name, "/cell_info_", velm_3rd_trim@project.name, ".csv"))
cell_info$X <- NULL


df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_3rd_trim_scenic)


