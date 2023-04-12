library(Seurat)
library(SeuratObject)
library(patchwork)
library(stringr)
library(dplyr)


main <- "/Home/ii/auraz/data/UCSC/outputs/Eze_Nowakowski_Velmeshev_2nd_trimester_integrated"

dir.create(main, recursive = T, showWarnings = F)

# Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html

input_rds_path <-  "/Home/ii/auraz/data/UCSC/Seurat_UCSC/integrated"
input_rds_files <- c("/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/Eze_2021.rds",
                     "/Home/ii/auraz/data/UCSC/Seurat_UCSC/others/Nowakowski_2017.rds",
                     "/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_2nd_trimester.rds")
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- input_rds_files
names(input_rds) <- str_remove_all(names(input_rds), ".rds")
names(input_rds) <- str_extract(names(input_rds), "\\w+$")

# normalize and identify variable features for each dataset independently
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Next, select features for downstream integration, and run PCA on each object in the list, which is required for running the alternative reciprocal PCA workflow
features <- SelectIntegrationFeatures(object.list = input_rds)
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# this command creates an 'integrated' data assay
anchors <- FindIntegrationAnchors(object.list = input_rds, reduction = "rpca", dims = 1:50)
#anchors <- FindIntegrationAnchors(object.list = input_rds, reference = c(1, 2), reduction = "rpca",dims = 1:50)
trim_2nd <- IntegrateData(anchorset = anchors, dims = 1:50)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(trim_2nd) <- "integrated"

# Run the standard workflow for visualization and clustering
trim_2nd <- ScaleData(trim_2nd, verbose = T)
trim_2nd <- RunPCA(trim_2nd, npcs = 30, verbose = T)
trim_2nd <- RunUMAP(trim_2nd, reduction = "pca", dims = 1:30)
trim_2nd <- FindNeighbors(trim_2nd, reduction = "pca", dims = 1:30)
trim_2nd <- FindClusters(trim_2nd, resolution = 0.5)

trim_2nd@project.name <- "trim_2nd_all"

trim_2nd@meta.data$cluster_final <- coalesce(trim_2nd@meta.data$Cell.Type, trim_2nd@meta.data$WGCNAcluster, trim_2nd@meta.data$cluster_final)
trim_2nd@meta.data[which(trim_2nd@meta.data$cluster_final== ""), "cluster_final"] <- "U2"

cluster_match <- c(
  "IPC-nEN1" = "IPC",     
  "IPC-div2" =   "IPC",   
  "IPC-div1" =  "IPC",     
  "OPC"     =  "Other",     
  "tRG"   =  "Radial Glial",       
  "IN-CTX-MGE2" = "Neuronal", 
  "nEN-late" = "Neuronal",    
  "nEN-early2" = "Neuronal", 
  "Mural"   =  "Other",          
  "IN-CTX-MGE1" = "Neuronal", 
  "IN-CTX-CGE2" = "Neuronal", 
  "RG-div2"  =  "Radial Glial",       
  "vRG"    =  "Radial Glial",         
  "IPC-nEN2"   = "IPC",       
  "IPC-nEN3"   = "IPC",       
  "IN-CTX-CGE1"= "Neuronal", 
  "IN-STR"     = "Neuronal",  
  "nIN5"    = "Neuronal",     
  "MGE-IPC3"     = "IPC",      
  "nIN4"       = "Neuronal",   
  "MGE-RG1"   =  "Radial Glial",     
  "MGE-RG2"  =  "Radial Glial",      
  "MGE-div"     =  "Other",    
  "MGE-IPC2"    = "IPC",   
  "MGE-IPC1"   = "IPC",     
  "nIN2"     = "Neuronal",      
  "nIN1"      = "Neuronal",     
  "EN-V1-1"     = "Neuronal",   
  "nIN3"       = "Neuronal",    
  "EN-PFC1"      = "Neuronal",  
  "nEN-early1"   = "Neuronal",  
  "U2"            = "Unknown", 
  "RG-div1"    =  "Radial Glial",  
  "EN-PFC2"     = "Neuronal", 
  "EN-PFC3"    = "Neuronal",  
  "EN-V1-2"    = "Neuronal",  
  "RG-early"    =  "Radial Glial",   
  "Endothelial"  =  "Other",     
  "Astrocyte"    =  "Other",     
  "Microglia"   =  "Other",     
  "oRG"         =  "Radial Glial",   
  "Glyc"        =  "Other",    
  "EN-V1-3"      = "Neuronal",  
  "Choroid"      =  "Other",     
  "U1"          = "Unknown", 
  "U3"          = "Unknown", 
  "U4"= "Unknown"
)

for (nowa_ct in names(cluster_match)) {
  trim_2nd@meta.data[which(trim_2nd@meta.data$cluster_final==nowa_ct), "cluster_final"] <- cluster_match[[nowa_ct]]
}

trim_2nd@meta.data$age_final <- round(coalesce(trim_2nd@meta.data$Age,trim_2nd@meta.data$Age_in_Weeks))
trim_2nd@meta.data[which(is.na(trim_2nd@meta.data$age_final)), "age_final"] <- "2nd"

trim_2nd@meta.data$trimester <- rep("trim", nrow(trim_2nd@meta.data))

trim_2nd@meta.data[which(0 < as.numeric(trim_2nd@meta.data$age_final) & as.numeric(trim_2nd@meta.data$age_final) < 13), "trimester"] <- "1st"
trim_2nd@meta.data[which(12 < as.numeric(trim_2nd@meta.data$age_final) & as.numeric(trim_2nd@meta.data$age_final) < 27), "trimester"] <- "2nd"
trim_2nd@meta.data[which(26 < as.numeric(trim_2nd@meta.data$age_final) & as.numeric(trim_2nd@meta.data$age_final) < 41), "trimester"] <- "3rd"

trim_2nd@meta.data[which(trim_2nd@meta.data$trimester=="trim"), "trimester"] <- trim_2nd@meta.data[which(trim_2nd@meta.data$trimester=="trim"), "age_final"]

# Visualization
pdf(paste0(main, "/", trim_2nd@project.name, "_trimesters.pdf"))
DimPlot(trim_2nd, reduction = "umap", group.by = "trimester")
dev.off()

pdf(paste0(main, "/", trim_2nd@project.name, "_cluster_final.pdf"))
DimPlot(trim_2nd, reduction = "umap", group.by = "cluster_final")
dev.off()

saveRDS(trim_2nd, paste0(input_rds_path, "/", trim_2nd@project.name, ".rds"))


########

trim_2nd <- readRDS(paste0(input_rds_path, "/trim_2nd_all.rds"))

pdf(paste0(main, "/", trim_2nd@project.name, "_projects.pdf"))
DimPlot(trim_2nd, reduction = "umap", group.by = "proj")
dev.off()

