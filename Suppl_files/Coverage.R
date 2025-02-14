# Author: Aura Zelco
# Brief description:
# This script is used for plotting the mapping coverage and gene coverage in each study
  # Brief procedure:
    # 1. imports the necessary libraries and defines a list of the Seurat Objects and the relative directories
    # 2. for each object, imports the RDS, creates the nFeatureRNA and nCountRNA plots and prints them as a pdf into the output path
  # Documentation abbreviations:
    # obj: object

#---------------------------------------------------------------------------------------------------
####### ON FURU

# 0. Import Libraries

library(Seurat)

main_furu <- "data/UCSC/Seurat_UCSC/Velmeshev/"
out_path <- "data/UCSC/Seurat_UCSC/Extra_figures/Gene_coverage/"
dir.create(out_path, recursive = T, showWarnings = F)

ds_paths <- c("Velmeshev_2nd_trimester"="Velmeshev_2022_2nd_trimester",
              "Velmeshev_10_20_years"="Velmeshev_2022_10_20_years")

ExtractMetadata <- function(main_dir, obj_path) {
  metadata_ls <- list()
  for (id in names(obj_path)) {
    print(id)
    obj_id <- readRDS(paste0(main_dir, obj_path[[id]], ".rds"))
    metadata_ls <- append(metadata_ls, list(obj_id@meta.data))
    rm(obj_id)
  }
  names(metadata_ls) <- names(obj_path)
  return(metadata_ls)
}

furu_meta_ls <- ExtractMetadata(main_furu, ds_paths)

saveRDS(furu_meta_ls, paste0(out_path, "furu_metadata_coverage.rds"))

# transfer to local machine using sftp

#---------------------------------------------------------------------------------------------------
####### LOCAL

# 0. Import Libraries

library(Seurat)
library(ggplot2)
library(ggpubr)

main <- "data/"


# 1. Read local RDS into a list

ds_paths <- c("Velmeshev_3rd_trimester"= "UCSC/Seurat_UCSC/",
              "Velmeshev_0_1_years"="UCSC/Seurat_UCSC/",
              "Velmeshev_1_2_years"= "UCSC/Seurat_UCSC/",
              "Velmeshev_2_4_years"= "UCSC/Seurat_UCSC/",
              "Velmeshev_Adults" = "UCSC/Seurat_UCSC/",
              "GSE157827_Healthy" = "DISCOv1.0/individual_proj/",
              "GSE157827_AD" = "DISCOv1.0/individual_proj/",
              "GSE174367_Healthy" = "DISCOv1.0/individual_proj/",
              "GSE174367_AD" = "DISCOv1.0/individual_proj/",
              "PRJNA544731_Healthy" = "DISCOv1.0/individual_proj/",
              "PRJNA544731_MS" = "DISCOv1.0/individual_proj/"
              )

ds_names <- c("Velmeshev_3rd_trimester"= "Velmeshev_2022_3rd_trimester",
              "Velmeshev_0_1_years"= "Velmeshev_2022_0_1_years",
              "Velmeshev_1_2_years"= "Velmeshev_2022_1_2_years",
              "Velmeshev_2_4_years"= "Velmeshev_2022_2_4_years",
              "Velmeshev_Adults"= "Velmeshev_2022_Adult",
              "GSE157827_Healthy" = "GSE157827_Healthy",
              "GSE157827_AD" = "GSE157827_AD",
              "GSE174367_Healthy" = "GSE174367_Healthy",
              "GSE174367_AD" = "GSE174367_AD",
              "PRJNA544731_Healthy" = "PRJNA544731_Healthy" ,
              "PRJNA544731_MS" = "PRJNA544731_MS" 
              )

out_path <- paste0(main, "Extra_figures/Gene_coverage/")
dir.create(out_path)


ExtractMetadata <- function(main_dir, obj_name, obj_path) {
  metadata_ls <- list()
  for (id in names(obj_path)) {
    print(id)
    obj_id <- readRDS(paste0(main_dir, obj_path[[id]], obj_name[[id]], ".rds"))
    metadata_ls <- append(metadata_ls, list(obj_id@meta.data))
    rm(obj_id)
  }
  names(metadata_ls) <- names(obj_path)
  return(metadata_ls)
}

ExtractParams <- function(metadata_ls, params, ds_order) {
  features_ls <- list()
  for (id in names(metadata_ls)) {
    features_ls <- append(features_ls, list(metadata_ls[[id]][, params]))
  }
  names(features_ls) <- names(metadata_ls)
  features_df <- do.call("rbind", features_ls)
  features_df <- cbind("ds" = gsub("\\..*","", rownames(features_df)), features_df)
  rownames(features_df) <- NULL
  features_df$ds <- factor(features_df$ds, ds_order)
  return(features_df)
}

meta_ls <- ExtractMetadata(main, ds_names, ds_paths)

furu_meta_ls <- readRDS(paste0(main, "Extra_figures/Gene_coverage/furu_metadata_coverage.rds"))

order_ds <- c("Velmeshev_2nd_trimester",
              "Velmeshev_3rd_trimester", 
              "Velmeshev_0_1_years", 
              "Velmeshev_1_2_years", 
              "Velmeshev_2_4_years",
              "Velmeshev_10_20_years",
              "Velmeshev_Adults",
              "GSE157827_Healthy",
              "GSE174367_Healthy",
              "PRJNA544731_Healthy",
              "GSE157827_AD",
              "GSE174367_AD",
              "PRJNA544731_MS")

feat_df <- ExtractParams(c(meta_ls, furu_meta_ls), c("nFeature_RNA", "nCount_RNA"), order_ds)

Feature_RNA <- function(features_df, legend_shown="bottom") {
  violin_plt <- ggplot(features_df, aes(x = ds, y = nFeature_RNA, fill = ds)) + 
      geom_violin() +
      #geom_jitter(shape=16, position=position_jitter(0.2), color="black", size=0.1, alpha = 0.2) +
      labs(fill = "Datasets", title="A)") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            strip.text.y = element_text(size=8, face="bold", colour = "black", angle = 0),
            strip.text.x = element_text(size=10, face="bold", colour = "black"),
            plot.title = element_text(size=12, face="bold", colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position = legend_shown, 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  return(violin_plt)
}

Count_RNA <- function(features_df) {
  violin_plt <- ggplot(features_df, aes(x = ds, y = nCount_RNA, fill = ds)) + 
    geom_violin() +
    #geom_jitter(shape=16, position=position_jitter(0.2), color="black", size=0.1, alpha = 0.2) +
    labs(fill = "Datasets", title="B)") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text.y = element_text(size=8, face="bold", colour = "black", angle = 0),
          strip.text.x = element_text(size=10, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(violin_plt)
}


FeatureScatterPlot <- function(features_df) {
  scatter_plt <- ggplot(features_df, aes(x = nCount_RNA, y = nFeature_RNA, color = ds)) + 
    geom_point() +
    labs(fill = "Datasets") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text.y = element_text(size=8, face="bold", colour = "black", angle = 0),
          strip.text.x = element_text(size=10, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(scatter_plt)
}


feat <- Feature_RNA(feat_df, legend_shown = "none")
count <- Count_RNA(feat_df)



pdf(paste0(out_path, "nFeature_RNA_ds.pdf"), width = 10, height = 8)
print(Feature_RNA(feat_df))
dev.off()

pdf(paste0(out_path, "nCount_RNA_ds.pdf"), width = 10, height = 8)
print(Count_RNA(feat_df))
dev.off()


pdf(paste0(out_path, "RNA_ds_features.pdf"), width = 10, height = 12)
print(ggarrange(feat, count, nrow=2))
dev.off()


