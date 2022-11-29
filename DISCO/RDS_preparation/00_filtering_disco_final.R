library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)
library(matrixStats)
library(stringr)
library(purrr)
`%!in%` <- Negate(`%in%`)

disco_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/"
maindir <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/"


##################################### Filter original DISCO v1.0 dataset

# load the dataset
disco_brain <- readRDS(paste0(disco_path, "brainV1.0.rds"))

# plot initial metadata
DimPlot(disco_brain, reduction = "umap")
DimPlot(disco_brain, reduction = "umap", group.by = "project_id")
DimPlot(disco_brain, reduction = "umap", group.by = "sample_type") # -> normal
DimPlot(disco_brain, reduction = "umap", group.by = "anatomical_site") # remove NA and midbrain
DimPlot(disco_brain, reduction = "umap", group.by = "disease")
DimPlot(disco_brain, reduction = "umap", group.by = "tissue") # all brain
DimPlot(disco_brain, reduction = "umap", group.by = "platform") # 10x3'v2 and 10x3'v3
DimPlot(disco_brain, reduction = "umap", group.by = "gender") # F, M, NA
DimPlot(disco_brain, reduction = "umap", group.by = "race") # Caucasian, NA

allmisscols <- sapply(disco_brain@meta.data, function(x) all(is.na(x) | x == '' ))
allmisscols # -> these columns have NAs, so cannot be plotted

FeaturePlot(disco_brain, features = "XIST", split.by = "gender")

disco_FM <- subset(disco_brain, subset = gender ==  c("F","M"))

DimPlot(disco_FM, reduction = "umap", group.by = "gender")
DimPlot(disco_FM, reduction = "umap", group.by = "age")


samples_cell_num <- as.data.frame(table(disco_FM$sample))
colnames(samples_cell_num) <- c("samples", "count")
samples_cell_num$percent <- samples_cell_num$count * 100 / sum(samples_cell_num$count)

samples_cell_num$sex <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$age <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$proj <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$platform <- rep(NA, length(samples_cell_num$samples))



for (sample in samples_cell_num$samples) {
  samples_cell_num[which(samples_cell_num$samples==sample), "proj"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "project_id"])
  samples_cell_num[which(samples_cell_num$samples==sample), "sex"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "gender"])
  samples_cell_num[which(samples_cell_num$samples==sample), "age"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "age"])
  samples_cell_num[which(samples_cell_num$samples==sample), "platform"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "platform"])
}

col_factors <- c("samples", 
                 "sex",                          
                 "age",
                 "proj"
)

samples_cell_num[col_factors] <- lapply(samples_cell_num[col_factors], as.factor)  

ggplot(samples_cell_num,aes(samples, count, fill=proj))+
  geom_bar(stat="identity") +
  facet_wrap(~sex*platform, scales="free_x", dir="v") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

sum(samples_cell_num$sex=="F") # 28
sum(samples_cell_num$sex=="M") # 75

ggsave(paste0(disco_path, "DISCO_FM/samples_counts.png"))

FeaturePlot(disco_FM, features = "XIST", split.by = "gender", reduction = "umap")

ggsave(paste0(disco_path, "DISCO_FM/xist_exp.png"))

celltype_sample <- paste(disco_FM@meta.data$sample, disco_FM@meta.data$ct, sep="_")

disco_FM@meta.data$ct_sample <- celltype_sample

celltype_num <- as.data.frame(table(disco_FM$ct_sample))
colnames(celltype_num) <- c("name", "freq")
library(tidyr)
celltype_num <- separate(celltype_num, name, into = c("sample", "ct"), sep = "_")


celltype_num$sex <- rep(NA, length(celltype_num$sample))
celltype_num$age <- rep(NA, length(celltype_num$sample))
celltype_num$proj <- rep(NA, length(celltype_num$sample))
celltype_num$platform <- rep(NA, length(celltype_num$sample))
celltype_num$anatomic_tissue <- rep(NA, length(celltype_num$sample))

for (sample in celltype_num$sample) {
  celltype_num[which(celltype_num$sample==sample), "proj"] <- unique(disco@meta.data[which(disco@meta.data$sample == sample), "project_id"])
  celltype_num[which(celltype_num$sample==sample), "sex"] <- unique(disco@meta.data[which(disco@meta.data$sample == sample), "gender"])
  celltype_num[which(celltype_num$sample==sample), "age"] <- unique(disco@meta.data[which(disco@meta.data$sample == sample), "age"])
  celltype_num[which(celltype_num$sample==sample), "platform"] <- unique(disco@meta.data[which(disco@meta.data$sample == sample), "platform"])
  celltype_num[which(celltype_num$sample==sample), "anatomic_tissue"] <- unique(disco@meta.data[which(disco@meta.data$sample == sample), "anatomical_site"])
}

col_factors <- c("sample", 
                 "ct",
                 "sex",                          
                 "age",
                 "proj",
                 "platform",
                 "anatomic_tissue"
)

celltype_num[col_factors] <- lapply(celltype_num[col_factors], as.factor)  

write.csv(celltype_num, paste0(disco_path, "celltype_num.csv"))

ggplot(celltype_num,aes(sample, freq, fill=proj))+
  geom_bar(stat="identity") +
  facet_wrap(~sex*platform, scales="free_x", dir="v") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

# same as before!

ggplot(celltype_num, aes(ct, freq, fill=ct)) + 
  geom_bar(stat="identity", position="stack") +
  facet_wrap(~sex, scales="free_x") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

ggsave(paste0(disco_path, "DISCO_FM/ct_counts.png"))

saveRDS(disco_FM, paste0(disco_path, "DISCOv1.0_FM/disco_FM.rds"))

disco <- readRDS(paste0(disco_path, "brainV1.0_all_FM.rds"))

celltype_num <- read.csv(paste0(disco_path, "celltype_num.csv"))
celltype_num[,1] <- NULL
col_factors <- c("sample", 
                 "ct",
                 "sex",                          
                 "age",
                 "proj",
                 "platform",
                 "anatomic_tissue"
)
celltype_num[col_factors] <- lapply(celltype_num[col_factors], as.factor)  


keep_proj <- vector()

for (proj in levels(celltype_num$proj)) {
  if (length(unique(celltype_num[which(celltype_num$proj==proj), "sex"])) == 2) {
    keep_proj <- c(keep_proj, proj)
  }
}

disco_filt <- subset(disco, subset =  project_id %in% keep_proj & anatomical_site != "midbrain")

rm(disco)

celltype_num <- as.data.frame(table(disco_filt$ct_sample))
colnames(celltype_num) <- c("name", "freq")
library(tidyr)
celltype_num <- separate(celltype_num, name, into = c("sample", "ct"), sep = "_")


celltype_num$sex <- rep(NA, length(celltype_num$sample))
celltype_num$age <- rep(NA, length(celltype_num$sample))
celltype_num$proj <- rep(NA, length(celltype_num$sample))
celltype_num$platform <- rep(NA, length(celltype_num$sample))
celltype_num$anatomic_tissue <- rep(NA, length(celltype_num$sample))

for (sample in celltype_num$sample) {
  celltype_num[which(celltype_num$sample==sample), "proj"] <- unique(disco_filt@meta.data[which(disco_filt@meta.data$sample == sample), "project_id"])
  celltype_num[which(celltype_num$sample==sample), "sex"] <- unique(disco_filt@meta.data[which(disco_filt@meta.data$sample == sample), "gender"])
  celltype_num[which(celltype_num$sample==sample), "age"] <- unique(disco_filt@meta.data[which(disco_filt@meta.data$sample == sample), "age"])
  celltype_num[which(celltype_num$sample==sample), "platform"] <- unique(disco_filt@meta.data[which(disco_filt@meta.data$sample == sample), "platform"])
  celltype_num[which(celltype_num$sample==sample), "anatomic_tissue"] <- unique(disco_filt@meta.data[which(disco_filt@meta.data$sample == sample), "anatomical_site"])
}

col_factors <- c("sample", 
                 "ct",
                 "sex",                          
                 "age",
                 "proj",
                 "platform",
                 "anatomic_tissue"
)

celltype_num[col_factors] <- lapply(celltype_num[col_factors], as.factor)  

write.csv(celltype_num, paste0(disco_path, "celltype_num_filt.csv"))

ggplot(celltype_num,aes(sample, freq, fill=proj))+
  geom_bar(stat="identity") +
  facet_wrap(~sex*platform, scales="free_x", dir="v") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

ggplot(celltype_num,aes(sample, freq, fill=anatomic_tissue))+
  geom_bar(stat="identity") +
  facet_wrap(~sex, scales="free_x", dir="v") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

project_sex <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, sep="_")
disco_filt@meta.data$proj_sex <- project_sex

proj_sex_ct <- paste(disco_filt@meta.data$proj_sex, disco_filt@meta.data$ct, sep="_")
disco_filt@meta.data$proj_sex_ct <- proj_sex_ct

saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

disco_filt@meta.data$sex_ct <- paste(disco_filt@meta.data$gender, disco_filt@meta.data$ct, sep="_")
disco_filt@meta.data[which(is.na(disco_filt@meta.data$disease)), "disease"] <- "Normal"
disco_filt@meta.data$proj_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data$proj_sex_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data[which(disco_filt@meta.data$disease=="Alzheimer’s disease"), "disease"] <- "Alzheimer's disease"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="multiple sclerosis"), "disease"] <- "Multiple Sclerosis"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="neurosurgical disease"), "disease"] <- "Neurosurgical Disease"  
disco_filt@meta.data$proj_sex_disease_ct <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, disco_filt@meta.data$ct, sep="_")


##################################### Filter cts for DEGs analysis steps

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

num_proj_sex_disease_ct <- as.data.frame(table(disco_filt$proj_sex_disease_ct))
num_proj_sex_disease_ct <- separate(num_proj_sex_disease_ct, Var1, into = c("proj", "sex" , "disease", "ct"), sep = "_")

col_factors <- c("proj", 
                 "sex",
                 "disease",
                 "ct"
)
num_proj_sex_disease_ct[col_factors] <- lapply(num_proj_sex_disease_ct[col_factors], as.factor)  
names(num_proj_sex_disease_ct)[names(num_proj_sex_disease_ct) == 'Freq'] <- "count"

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

num_normal <- subset(num_proj_sex_disease_ct, disease == "Normal")
num_AD <- subset(num_proj_sex_disease_ct, disease == "Alzheimer's disease")
num_MS <- subset(num_proj_sex_disease_ct, disease == "Multiple Sclerosis")

min_num_cells <- c(10,50,100)

for (min_cells in min_num_cells) {
  num_normal_filt <- FiltDF(num_normal, "Normal", min_cells)
  num_AD_filt <- FiltDF(num_AD, "AD", min_cells)
  num_MS_filt <- FiltDF(num_MS, "MS", min_cells)
  num_filt <- rbind(num_normal_filt, num_AD_filt, num_MS_filt)
  num_filt$idents <- paste(num_filt$proj, num_filt$sex, num_filt$disease, num_filt$ct, sep="_")
  num_filt$name_subfolders <- str_replace_all(num_filt$ct, "/", "_")
  #print(nrow(num_filt))
  write.csv(num_filt, file = paste0(disco_path, 
                                        "DEGs/outputs/final_filt_", min_cells, ".csv"),
            row.names = F)
}

for (min_cells in min_num_cells) {
  for (id in levels(num_proj_sex_disease_ct$proj)) {
    pdf(paste0(disco_path, "DEGs/outputs/", id, "_filt_counts_", min_cells, ".pdf"), 10, 15)
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

saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

##################################### INFO FOR SCENIC PIPELINE

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

num_proj_sex <- as.data.frame(table(disco_filt$proj_sex_disease))
num_proj_sex <- separate(num_proj_sex, Var1, into = c("proj", "sex", "disease"), sep = "_")
names(num_proj_sex)[names(num_proj_sex) == 'Freq'] <- "count"
write.csv(num_proj_sex, paste0(disco_path, "DEGs/num_proj_sex_disease.csv"))

for (dis_type in unique(num_proj_sex$disease)) {
  pdf(paste0(disco_path, "DEGs/num_proj_sex_", dis_type,".pdf"))
  print(
    ggplot(num_proj_sex[which(num_proj_sex$disease==dis_type), ], aes(proj, count, fill=sex)) +
      geom_bar(stat="identity", position="dodge") +
      labs(x="Project ID", y='Cell Counts', fill="Sex", main = dis_type) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.5, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom")
  )
  dev.off()
}

num_proj_sex_ct <- as.data.frame(table(disco_filt$proj_sex_disease_ct))
num_proj_sex_ct <- separate(num_proj_sex_ct, Var1, into = c("proj", "sex", "disease", "ct"), sep = "_")
names(num_proj_sex_ct)[names(num_proj_sex_ct) == 'Freq'] <- "count"
write.csv(num_proj_sex_ct, paste0(disco_path, "DEGs/num_proj_sex_ct.csv"))

for (dis_type in unique(num_proj_sex_ct$disease)) {
  pdf(paste0(disco_path, "DEGs/num_proj_sex_", dis_type,"_ct.pdf"),
      width = 15, height = 18.75)
  print(
    ggplot(num_proj_sex_ct[which(num_proj_sex_ct$disease==dis_type), ], aes(proj, count, fill=sex)) +
      geom_bar(stat="identity", position="dodge") +
      labs(x="Project ID", y='Cell Counts', fill="Sex", main = dis_type) +
      facet_wrap(~ct, scales = "free") +
      geom_hline(yintercept = 100, linetype=2) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.5, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom")
  )
  dev.off()
}

saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

######## HIGHEST VARIABLE GENES - 2022.10.18

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

disco_filt@meta.data$proj_sex_disease_ct_sample <- paste(disco_filt@meta.data$proj_sex_disease_ct, disco_filt@meta.data$sample_id, sep="_")

Idents(disco_filt) <- "proj_sex_disease_ct_sample"

expr_mat_all <- GetAssayData(disco_filt[["RNA"]], slot="data")

Idents(disco_filt) <- "proj_sex_disease_ct"

cell_info <- data.frame()

for (i in unique(disco_filt@meta.data$proj_sex_disease_ct)) {
  print(i)
  cell_id <- WhichCells(disco_filt, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}

cell_info <- separate(cell_info, cell_id, into = c("barcode", "sample"), sep = "--", remove = FALSE)

write.csv(cell_info,
          paste0(disco_path, "DEGs/cell_info.csv"))

saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

rm(disco_filt)

expr_mat_all <- as.data.frame(as.matrix(expr_mat_all))
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

saveRDS(expr_mat_all, paste0(disco_path, "DEGs/top_2000_SD_expr_matrix.rds"))

##### Map the samples back to the groups they belong to
cell_info$cell_id <- paste(cell_info$barcode, cell_info$sample, sep = "--")
cell_info <- cell_info %>% 
  relocate(cell_id, .before = barcode)

expr_mat_all <- readRDS(paste0(disco_path, "2SCENIC/top_2000_SD_expr_matrix.rds"))

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
    if (grepl("F", i)) {
      m_id <- str_replace(i, "F", "M")
      if (m_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, m_id)
      }
    } else {
      f_id <- str_replace(i, "M", "F")
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

remove_groups <- vector()
if (length(group_list) %% 2 != 0 ) {
  f_list <- vector()
  m_list <- vector()
  for (i in names(group_list)) {
    if (grepl("F", i)) {
      gen <- str_remove(i, "F_")
      f_list <- c(f_list, gen)
    } else {
      gen <- str_remove(i, "M_")
      m_list <- c(m_list, gen)
    }
  }
  if (identical(m_list, f_list) == FALSE) {
    if (length(m_list) > length(f_list)) {
      remove_groups <- c(m_list[which(m_list %!in% f_list)], "M")
    } else {
      remove_groups <- c(f_list[which(f_list %!in% m_list)], "F")
    }
  }
}

# "PRJNA544731_M_Normal_microglia"
group_list[["PRJNA544731_M_Normal_microglia"]] <- NULL

for (id in names(group_list)){
  if (ncol(group_list[[id]]) > 101) {
    print(ncol(group_list[[id]]))
    print(ncol(group_list[[id]]) %/% 100 )
  }
}


ncol(group_list[[1]])
for (k in 1:(ncol(group_list[[1]]) %/% 100)) {
  print(k)
}

hist(sapply(1:length(names(group_list)), function(i) ncol(group_list[[i]])),
     breaks = 150,
     xlim=c(0,1000))


group_list100 <- group_list
group_list500 <- group_list

group_list100 <- remove_dfs(group_list100, 100)
group_list500 <- remove_dfs(group_list500, 500)


plot_group_numbers <- function(df_list, thresh) {
  ids <- as.data.frame(names(df_list))
  colnames(ids) <- c("groups")
  ids <- separate(ids, groups, into = c("proj", "sex", "disease", "ct"), sep="_", remove=FALSE)
  col_factors <- c("proj", "sex", "disease", "ct")
  ids[col_factors] <- lapply(ids[col_factors], as.factor)  
  ids$length_groups <- sapply(1:length(names(df_list)), function(i) ncol(df_list[[i]]))
  pdf(paste0(disco_path, "SCENIC/num_filt_", thresh, "_cells.pdf"),
      height = 12)
  print(
    ggplot(ids, aes(ct, length_groups, fill=sex)) +
      geom_bar(stat="identity", position = "dodge") + 
      facet_wrap(~proj*disease, nrow = 3, scales = "free") +
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

plot_group_numbers(group_list, 10)
plot_group_numbers(group_list100, 100)
plot_group_numbers(group_list500, 500)

###### Create Randomly sampled dfs

expr_mat_all <- readRDS(paste0(disco_path, "SCENIC/top_2000_SD_expr_matrix.rds"))
cell_info <- read.csv(paste0(disco_path, "SCENIC/cell_info.csv"))

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
      m_id <- str_replace(i, "F", "M")
      if (m_id %!in% incomplete_dfs) {
        add_counterpart <- c(add_counterpart, m_id)
      }
    } else {
      f_id <- str_replace(i, "M", "F")
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

sub_disease <- list.dirs(maindir, recursive=FALSE, full.names = FALSE)

dfs_split <- list("Alzheimer's disease" = vector(), 
                  "Multiple Sclerosis" = vector(),
                  "Normal"   = vector())

for (dis_type in sub_disease) {
  dfs_split[[dis_type]] <- names(df_list100)[which(grepl(dis_type, names(df_list100)))]
}


norm <- df_list100[dfs_split[["Normal"]]]
ad <- df_list100[dfs_split[["Alzheimer's disease"]]]
ms <- df_list100[dfs_split[["Multiple Sclerosis"]]]

check_dfs <- function(group_list) {
  remove_groups <- vector()
  if (length(group_list) %% 2 != 0 ) {
    f_list <- vector()
    m_list <- vector()
    for (i in names(group_list)) {
      if (grepl("F", i)) {
        gen <- str_remove(i, "F_")
        f_list <- c(f_list, gen)
      } else {
        gen <- str_remove(i, "M_")
        m_list <- c(m_list, gen)
      }
    }
    if (identical(m_list, f_list) == FALSE) {
      if (length(m_list) > length(f_list)) {
        remove_groups <- c(m_list[which(m_list %!in% f_list)], "M")
      } else {
        remove_groups <- c(f_list[which(f_list %!in% m_list)], "F")
      }
    }
  }
  return(remove_groups)
}

check_dfs(norm)
check_dfs(ad)
check_dfs(ms)

rand_sample <- function(group_list, num_sampling, num_cells, main, dis_type) {
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
  names(sampled_dfs) <- lapply(1:length(sampled_names), function(i) str_replace_all(sampled_names[i],
                                                                                    "/", "_"))
  names(sampled_dfs) <- lapply(1:length(names(sampled_dfs)), function(i) str_replace_all(names(sampled_dfs)[i],
                                                                                         " ", "_"))
  dir.create(paste0(main, dis_type, "/sampled_", num_cells, "_cells"), showWarnings = FALSE)
  lapply(1:length(names(sampled_dfs)), function(i) write.csv(sampled_dfs[[i]], 
                                                             file = paste0(main, dis_type, "/sampled_", num_cells, "_cells/", names(sampled_dfs)[i], ".csv"),
                                                             row.names = FALSE))
  return(sampled_dfs)
}


norm_sampled <- rand_sample(norm, 3, 100, maindir, sub_disease[3])
ad_sampled <- rand_sample(ad, 3, 100, maindir, sub_disease[1])
ms_sampled <- rand_sample(ms, 3, 100, maindir, sub_disease[2])


############ For 02C_Conservation

disco_filt <- readRDS(paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

Idents(disco_filt) <- "proj_sex_disease_ct"

expr_mat_all_cts <- GetAssayData(disco_filt[["RNA"]], slot="data")
saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))

rm(disco_filt)
expr_mat_all_cts <- as.data.frame(as.matrix(expr_mat_all_cts))

cell_info <- read.csv(paste0(disco_path, "SCENIC/cell_info.csv"))
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

dfs_split <- list("Alzheimer'sdisease" = vector(), 
                  "Multiple Sclerosis" = vector(),
                  "Normal"   = vector())

sub_disease <- list.dirs(paste0(disco_path, "DEGs/outputs/"), recursive=FALSE, full.names = FALSE)

for (dis_type in sub_disease) {
  dfs_split[[dis_type]] <- names(df_list)[which(grepl(dis_type, names(df_list)))]
}


norm <- df_list[dfs_split[["Normal"]]]
ad <- df_list[dfs_split[["Alzheimer's disease"]]]
ms <- df_list[dfs_split[["Multiple Sclerosis"]]]


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

norm_df <- FiltDisDf(norm)
ad_df <- FiltDisDf(ad)
ms_df <- FiltDisDf(ms)

tot_df <- rbind(norm_df, ad_df, ms_df)

tot_df <- separate(tot_df, cts, into=c("proj", "sex", "disease", "ct"), sep ="_", remove = FALSE)
colnames(tot_df)
names(tot_df)[names(tot_df) == "cts"] <- "og"

col_factors <- c("og", "proj", "sex", "disease", "ct")
tot_df[col_factors] <- lapply(tot_df[col_factors], as.factor) 

diseases <- vector()
sexes <- vector()
cts <- vector()
genes <- vector()
for (dis_type in levels(tot_df$disease)) {
  if (dis_type == "Multiple Sclerosis") {
    for (sex_id in levels(tot_df$sex)) {
      for (ct_id in levels(tot_df$ct)) {
        proj_MS <- unique(tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id), "proj"])
        common_genes <- tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_MS), "genes"]
        genes <- c(genes, common_genes)
        diseases <- c(diseases, rep(dis_type, length(common_genes)))
        sexes <- c(sexes, rep(sex_id, length(common_genes)))
        cts <- c(cts, rep(ct_id, length(common_genes)))
      }
    }
  } else {
    for (sex_id in levels(tot_df$sex)) {
      for (ct_id in levels(tot_df$ct)) {
        proj_list <- unique(tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id), "proj"])
        if (length(proj_list)==3) {
          common_genes <- intersect(tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_list[1]), "genes"],
                                    intersect(
                                      tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_list[2]), "genes"],
                                      tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_list[3]), "genes"]
                                    ))
          genes <- c(genes, common_genes)
          diseases <- c(diseases, rep(dis_type, length(common_genes)))
          sexes <- c(sexes, rep(sex_id, length(common_genes)))
          cts <- c(cts, rep(ct_id, length(common_genes)))
        } else if (length(proj_list)==2) {
          common_genes <- intersect(tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_list[1]), "genes"],
                                    tot_df[which(tot_df$sex==sex_id & tot_df$disease==dis_type & tot_df$ct==ct_id & tot_df$proj==proj_list[2]), "genes"])
          genes <- c(genes, common_genes)
          diseases <- c(diseases, rep(dis_type, length(common_genes)))
          sexes <- c(sexes, rep(sex_id, length(common_genes)))
          cts <- c(cts, rep(ct_id, length(common_genes)))
        }
      }
    }
  }
}

tot_genes <- as.data.frame(cbind(diseases, sexes, cts, genes))
colnames(tot_genes) <- c("disease", "sex", "ct", "genes")
col_factors <- c("disease", "sex", "ct")
tot_genes[col_factors] <- lapply(tot_genes[col_factors], as.factor) 

write.csv(tot_genes, paste0(disco_path, "DEGs/tot_genes_ct.csv"))

