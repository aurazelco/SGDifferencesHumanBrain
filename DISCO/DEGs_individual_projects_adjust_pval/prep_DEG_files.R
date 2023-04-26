# Author: Aura Zelco
# Brief description:
  # This script is used to prepare files for the DEGs analysis
# Brief procedure: 
    # 1. filter cell types with less than 100 cells and plot cell number per sex, cell type and disease condition
    # 2. Save cell info (barcodes and associated metadata) and genes mtx for Conservation analysis

# OBS: this script requires visual inspection and manually-entered information, therefore should be opened in a R environment/IDE (e.g. RStudio). 

# Imports necessary libraries
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)

# Defines a custom function, used to subset using exclusion rather than inclusion criteria
`%!in%` <- Negate(`%in%`)

# Defines the main directories
disco_path <- ""

##################################### 1. Filter cts for DEGs analysis steps

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

# extract the number of cells per combination of project, sex, disease condition and celltype
num_proj_sex_disease_ct <- as.data.frame(table(disco_filt$proj_sex_disease_ct))
num_proj_sex_disease_ct <- separate(num_proj_sex_disease_ct, Var1, into = c("proj", "sex" , "disease", "ct"), sep = "_")

col_factors <- c("proj", 
                 "sex",
                 "disease",
                 "ct"
)
num_proj_sex_disease_ct[col_factors] <- lapply(num_proj_sex_disease_ct[col_factors], as.factor)  
names(num_proj_sex_disease_ct)[names(num_proj_sex_disease_ct) == 'Freq'] <- "count"

# filters the above df to check which celltypes contain less than a custom threshold number of cells
# then removes also the counterpart sex if only one of the 2 is less than the threshold
FiltDF <- function(df, disease, min_num_cells) {
  `%!in%` <- Negate(`%in%`)
  df <- droplevels(df)
  incomplete_proj <- vector()
  for (type in levels(df$ct)) {
    for (id in levels(df$proj)) {
      if ((nrow(subset(df, subset = (ct==type & proj==id)))%%2!=0) | (any(subset(df, subset = (ct==type & proj==id))[,5] < min_num_cells))) {
        incomplete_proj <- c(incomplete_proj, (paste(id, type, sep="_")))
      }
    }
  }
  df$og <- paste(df$proj, df$ct, sep="_")
  df_filt <- df[df$og %!in% incomplete_proj,]
  return(df_filt)
}

# splits the dfs in three subsets, one per disease condition
num_normal <- subset(num_proj_sex_disease_ct, disease == "Normal")
num_AD <- subset(num_proj_sex_disease_ct, disease == "Alzheimer's disease")
num_MS <- subset(num_proj_sex_disease_ct, disease == "Multiple Sclerosis")

min_num_cells <- c(10,50,100)

# saves the CSV outputs for the three thresholds
for (min_cells in min_num_cells) {
  num_normal_filt <- FiltDF(num_normal, "Normal", min_cells)
  num_AD_filt <- FiltDF(num_AD, "AD", min_cells)
  num_MS_filt <- FiltDF(num_MS, "MS", min_cells)
  num_filt <- rbind(num_normal_filt, num_AD_filt, num_MS_filt)
  num_filt$idents <- paste(num_filt$proj, num_filt$sex, num_filt$disease, num_filt$ct, sep="_")
  num_filt$name_subfolders <- str_replace_all(num_filt$ct, "/", "_")
  #print(nrow(num_filt))
  write.csv(num_filt, file = paste0(disco_path, 
                                    "DEGs_proj_adjust_pval/outputs/final_filt_", min_cells, ".csv"),
            row.names = F)
}

# and the corresponding plots
for (min_cells in min_num_cells) {
  for (id in levels(num_proj_sex_disease_ct$proj)) {
    pdf(paste0(disco_path, "DEGs_proj_adjust_pval/outputs/", id, "_filt_counts_", min_cells, ".pdf"), 10, 15)
    print(ggplot(num_proj_sex_disease_ct[which(num_proj_sex_disease_ct$proj==id),], aes(disease, count, fill=sex)) +
            geom_bar(stat="identity", position = "dodge") + 
            labs(x="", y="Nuclei count", fill="Sex") +
            facet_wrap(~ct, scales = "free") +
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
}

# updates the RDS
saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))


######## Cell info = barcodes and asscoiated metadata

# Extract the expression matrix to then later use in the SCENIC pipeline

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

disco_filt@meta.data$proj_sex_disease_ct_sample <- paste(disco_filt@meta.data$proj_sex_disease_ct, disco_filt@meta.data$sample_id, sep="_")

cell_info <- data.frame()

# creates a df containing metadata associated with the cell barcodes -> used later to create the randomly sampled matrices
for (i in unique(disco_filt@meta.data$proj_sex_disease_ct)) {
  print(i)
  cell_id <- WhichCells(disco_filt, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}

cell_info <- separate(cell_info, cell_id, into = c("barcode", "sample"), sep = "--", remove = FALSE)

write.csv(cell_info,
          paste0(disco_path, "DEGs_proj_adjust_pval/cell_info.csv"))

saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

rm(disco_filt)

##################################### 2. Save cell info (barcodes and associated metadata) and genes mtx for Conservation analysis

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

Idents(disco_filt) <- "proj_sex_disease_ct"

expr_mat_all_cts <- GetAssayData(disco_filt[["RNA"]], slot="data")
saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

rm(disco_filt)
expr_mat_all_cts <- as.data.frame(as.matrix(expr_mat_all_cts))

cell_info <- read.csv(paste0(disco_path, "DEGs_proj_adjust_pval/cell_info.csv"))
cell_info$X <- NULL

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all_cts[ , (names(expr_mat_all_cts) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

sub_projs <- list.dirs(paste0(disco_path, "DEGs_proj_adjust_pval/"), recursive=FALSE, full.names = FALSE)[-1]

dfs_split <- list()
for (proj_type in sub_projs) {
  dfs_split <- append(dfs_split, list(names(df_list)[which(grepl(proj_type, names(df_list)))]))
}
names(dfs_split) <- sub_projs


dfs_proj_dis <- list()
for (proj_id in sub_projs) {
  sub_disease <- list.dirs(paste0(disco_path, "DEGs_proj_adjust_pval/", proj_id), recursive=FALSE, full.names = FALSE)
  proj_list <- list()
  proj_tot <- df_list[dfs_split[[proj_id]]]
  for (dis_type in sub_disease) {
    proj_list <- append(proj_list, list(proj_tot[names(proj_tot)[which(grepl(dis_type, names(proj_tot)))]]))
  }
  names(proj_list) <- sub_disease
  dfs_proj_dis <- append(dfs_proj_dis, list(proj_list))
}
names(dfs_proj_dis) <- sub_projs

FiltDisDf <- function(df_list_dis) {
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

tot_df <- list()
for (proj_id in names(dfs_proj_dis)) {
  for (dis_type in names(dfs_proj_dis[[proj_id]])) {
    tot_df <- append(tot_df, list(FiltDisDf(dfs_proj_dis[[proj_id]][[dis_type]])))
  }
}
tot_df <- do.call(rbind, tot_df)

tot_df <- separate(tot_df, cts, into=c("proj", "sex", "disease", "ct"), sep ="_", remove = FALSE)
colnames(tot_df)
names(tot_df)[names(tot_df) == "cts"] <- "og"


write.csv(tot_df, paste0(disco_path, "DEGs_proj_adjust_pval/tot_genes_ct.csv"))


