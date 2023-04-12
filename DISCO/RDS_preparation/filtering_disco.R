# Author: Aura Zelco
# Brief description:
    # This script is used to filter the initial DISCO SeuratObject v1.0, as downloaded from the DISCO database
# Brief procedure: Explore the metadata initially present in the DISCO object,and filters the initial SeuratObject until only projects with cortical samples and both sexes are present

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
disco_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/"

##################################### 1. Filter original DISCO v1.0 dataset

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
DimPlot(disco_brain, reduction = "umap", group.by = "gender") # F, M, NA -> the database refers to the sex as gender -> if not in the original metadata, form hereafter we used the word sex
DimPlot(disco_brain, reduction = "umap", group.by = "race") # Caucasian, NA

allmisscols <- sapply(disco_brain@meta.data, function(x) all(is.na(x) | x == '' ))
allmisscols # -> these columns have NAs, so cannot be plotted

# explore XISt expression -> XIST is a gene expressed only in females
FeaturePlot(disco_brain, features = "XIST", split.by = "gender")

# subsets the object so only samples with known sex data are kept
disco_FM <- subset(disco_brain, subset = gender ==  c("F","M"))

# Plots UMAPs to explore the sex and the age groups
DimPlot(disco_FM, reduction = "umap", group.by = "gender")
DimPlot(disco_FM, reduction = "umap", group.by = "age")

# saves to a df the information about the cells per sample
samples_cell_num <- as.data.frame(table(disco_FM$sample))
colnames(samples_cell_num) <- c("samples", "count")
samples_cell_num$percent <- samples_cell_num$count * 100 / sum(samples_cell_num$count)

# creates NA columns which are then later replaced with the corresponding information
samples_cell_num$sex <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$age <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$proj <- rep(NA, length(samples_cell_num$samples))
samples_cell_num$platform <- rep(NA, length(samples_cell_num$samples))

# replacement of NAs with metadata
for (sample in samples_cell_num$samples) {
  samples_cell_num[which(samples_cell_num$samples==sample), "proj"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "project_id"])
  samples_cell_num[which(samples_cell_num$samples==sample), "sex"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "gender"])
  samples_cell_num[which(samples_cell_num$samples==sample), "age"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "age"])
  samples_cell_num[which(samples_cell_num$samples==sample), "platform"] <- unique(disco_FM@meta.data[which(disco_FM@meta.data$sample == sample), "platform"])
}

# changes column types to factor
col_factors <- c("samples", 
                 "sex",                          
                 "age",
                 "proj"
)
samples_cell_num[col_factors] <- lapply(samples_cell_num[col_factors], as.factor)  

# information about sequencing platform used and saves plot to the original directory
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
ggsave(paste0(disco_path, "DISCO_FM/samples_counts.png"))

# number of F and M in all projects
sum(samples_cell_num$sex=="F") # 28
sum(samples_cell_num$sex=="M") # 75

# generate FeaturePlot to confirm difference in expression of XIST 
FeaturePlot(disco_FM, features = "XIST", split.by = "gender", reduction = "umap")
ggsave(paste0(disco_path, "DISCO_FM/xist_exp.png"))

# creates a new metadata column
disco_FM@meta.data$ct_sample <- paste(disco_FM@meta.data$sample, disco_FM@meta.data$ct, sep="_")


# creates a new df, with information about how many cells per celltype 
celltype_num <- as.data.frame(table(disco_FM$ct_sample))
colnames(celltype_num) <- c("name", "freq")
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
# saves the df
write.csv(celltype_num, paste0(disco_path, "celltype_num.csv"))

# generates plot and saves it 
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

# saves the filtered object
saveRDS(disco_FM, paste0(disco_path, "DISCOv1.0_FM/disco_FM.rds"))


# imports filtered object
disco <- readRDS(paste0(disco_path, "brainV1.0_all_FM.rds"))

# imports the last df with the metadata information
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

# initializes an empty vector to contain the complete projects
keep_proj <- vector()

# for loop which checks which projects have both sexes, and saves the project id to the keep_project vector
for (proj in levels(celltype_num$proj)) {
  if (length(unique(celltype_num[which(celltype_num$proj==proj), "sex"])) == 2) {
    keep_proj <- c(keep_proj, proj)
  }
}

# subsets the Seurat Object again, keeping only the projects with both sexes and with sampling not from the midbrain -> only cortical samples are kept
disco_filt <- subset(disco, subset =  project_id %in% keep_proj & anatomical_site != "midbrain")

# removes the previous Seurat object to free memory
rm(disco)

# repeats the same as above, creates df to check number of cells per sample and per celltype
celltype_num <- as.data.frame(table(disco_filt$ct_sample))
colnames(celltype_num) <- c("name", "freq")
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


# creates new metadata, so the samples can be more easily divided
project_sex <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, sep="_")
disco_filt@meta.data$proj_sex <- project_sex

proj_sex_ct <- paste(disco_filt@meta.data$proj_sex, disco_filt@meta.data$ct, sep="_")
disco_filt@meta.data$proj_sex_ct <- proj_sex_ct

# fixes some discrepancies in the metadata
disco_filt@meta.data$sex_ct <- paste(disco_filt@meta.data$gender, disco_filt@meta.data$ct, sep="_")
disco_filt@meta.data[which(is.na(disco_filt@meta.data$disease)), "disease"] <- "Normal"
disco_filt@meta.data$proj_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data$proj_sex_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data[which(disco_filt@meta.data$disease=="Alzheimer’s disease"), "disease"] <- "Alzheimer's disease"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="multiple sclerosis"), "disease"] <- "Multiple Sclerosis"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="neurosurgical disease"), "disease"] <- "Neurosurgical Disease"  
disco_filt@meta.data$proj_sex_disease_ct <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, disco_filt@meta.data$ct, sep="_")

# saves the newly filtered object to anew RDS -> this is the final object we worked with from now on
saveRDS(disco_filt, paste0(disco_path, "brainV1.0_all_FM_filt.rds"))


