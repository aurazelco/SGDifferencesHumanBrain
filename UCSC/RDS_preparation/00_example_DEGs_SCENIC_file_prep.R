################ FOR DEGs ANALYSIS

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/UCSC/RDS_preparation/00_DEGs_and_SCENIC_files_prep.R")

velm_10_20_years <- readRDS(paste0(rds_path, "Velmeshev_2022_10_20_years.rds"))

velm_10_20_years@meta.data$sex_ct <- paste(velm_10_20_years@meta.data$sex, velm_10_20_years@meta.data$cluster_final, sep="_")

num_sex_ct <- NumSexCt(velm_10_20_years)

min_num_cells <- c(10,50,100)

velm_10_20_years_deg <- paste0(main_deg, velm_10_20_years@project.name, "/outputs/")
dir.create(velm_10_20_years_deg, recursive = T, showWarnings = F)

PlotFiltDf(min_num_cells, num_sex_ct, velm_10_20_years_deg)

saveRDS(velm_10_20_years, paste0(rds_path, velm_10_20_years@project.name, ".rds"))

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

df_list100 <- group_list
df_list100 <- RemoveDfs(df_list100, 100)

CheckDfs(df_list100)

RandomSampling(df_list100, 3, 100, velm_10_20_years_scenic)

saveRDS(velm_10_20_years, paste0(rds_path, velm_10_20_years@project.name, ".rds"))