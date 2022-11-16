library(Seurat)
library(SeuratData)
library(patchwork)
library(metap)
library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
library(scales)
library(readxl)
library(matrixStats)


main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/"


# Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html
input_rds_path <-  "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC"
input_rds_files <- list.files(path = input_rds_path, pattern = ".rds", full.names = T)[1:2]
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- list.files(path = input_rds_path, pattern = ".rds", full.names = F)[1:2]
names(input_rds) <- str_remove_all(names(input_rds), ".rds")

# normalize and identify variable features for each dataset independently
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = input_rds)
common.anchors <- FindIntegrationAnchors(object.list = input_rds, anchor.features = features)

# this command creates an 'integrated' data assay
rds.combined <- IntegrateData(anchorset = common.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(rds.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
rds.combined <- ScaleData(rds.combined, verbose = FALSE)
rds.combined <- RunPCA(rds.combined, npcs = 30, verbose = FALSE)
rds.combined <- RunUMAP(rds.combined, reduction = "pca", dims = 1:30)
rds.combined <- FindNeighbors(rds.combined, reduction = "pca", dims = 1:30)
rds.combined <- FindClusters(rds.combined, resolution = 0.5)

rds.combined@project.name <- "Eze_Nowakowski_integrated"

# Visualization
DimPlot(rds.combined, reduction = "umap")
DimPlot(rds.combined, reduction = "umap", group.by = "proj")
ggsave(paste0(main, rds.combined@project.name, "_proj.pdf"))

saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")
#rm(input_rds, input_rds_files, input_rds_path)


####### FIX AGE METADATA

DimPlot(rds.combined, reduction = "umap", group.by = "Age") # Eze
DimPlot(rds.combined, reduction = "umap", group.by = "Carnegie_Stage") # Eze
DimPlot(rds.combined, reduction = "umap", group.by = "Age.Range") # Eze
DimPlot(rds.combined, reduction = "umap", group.by = "sex_age") # shared
DimPlot(rds.combined, reduction = "umap", group.by = "Age_in_Weeks") # Nowakowksi
DimPlot(rds.combined, reduction = "umap", group.by = "Sample.Age..pcw.") # Nowakowksi

rds.combined@meta.data[which(rds.combined@meta.data$Age.Range=="Midlle"),"Age.Range"] <- "Middle"

sex_age_shared <- as.data.frame(table(rds.combined$sex_age))
sex_age_shared <- separate(sex_age_shared, Var1, into=c("sex", "age"), remove = F, sep = "_")
sex_age_shared$age <- as.numeric(sex_age_shared$age)
sex_age_shared$age_round <- round(sex_age_shared$age)

ggplot(sex_age_shared, aes(round(age))) +
  geom_histogram(binwidth = 1)

rds.combined@meta.data$age_final <- round(coalesce(rds.combined@meta.data$Age,rds.combined@meta.data$Age_in_Weeks))

DimPlot(rds.combined, reduction = "umap", group.by = "age_final")
ggsave(paste0(main, rds.combined@project.name, "_age_combined.pdf"))

# 1st_trimester: 1 - 12 
# 2nd_trimester: 13 - 26
# 3rd_trimester: 27 - 40

rds.combined@meta.data$trimesters <- rep("no_data", nrow(rds.combined@meta.data))

rds.combined@meta.data[which(0 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 13), "trimesters"] <- "1st"
rds.combined@meta.data[which(12 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 27), "trimesters"] <- "2nd"
rds.combined@meta.data[which(26 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 41), "trimesters"] <- "3rd"

DimPlot(rds.combined, reduction = "umap", group.by = "trimesters")
ggsave(paste0(main, rds.combined@project.name, "_trimesters.pdf"))


saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")

rds.combined@meta.data$sex_age <- paste(rds.combined@meta.data$sex, rds.combined@meta.data$age_final, sep="_")

DimPlot(rds.combined, reduction = "umap", group.by = "sex_age")

rds.combined@meta.data$id <- coalesce(rds.combined@meta.data$Individual, rds.combined@meta.data$Name)
rds.combined@meta.data$id_sex_age <- paste(rds.combined@meta.data$id, rds.combined@meta.data$sex_age, sep="_")

####### NUM OF CELLS BY AGE AND TRIMESTER

num_cells <- as.data.frame(table(rds.combined$id_sex_age))
num_cells <- separate(num_cells, Var1, into=c("id", "sex", "age"), sep="_")
num_cells$age <- as.numeric(num_cells$age) 

#trimesters <- data.frame(c(1,12.5,26.5), c(12.5,26.5,40.5), c(hue_pal()(3)))
#colnames(trimesters) <- c("bg_start", "bg_end", "colors")

p1 <- ggplot() +
  geom_bar(data = num_cells, aes(age, fill=sex)) +
  geom_vline(xintercept = c(12.5, 26.5), linetype = "dashed") +
  #geom_rect(data = trimesters, aes(xmin = bg_start,
  #          xmax = bg_end,
  #          ymin = - Inf,
  #          ymax = Inf,
  #          fill=colors),
  #          alpha = 0.1) +
  scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  scale_x_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="GW Age", y="Number of samples", fill="Sex") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

p2 <- ggplot(num_cells, aes(age, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  scale_x_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="GW Age", y="Number of cells/sample", fill="Sex") +
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

num_cells_int <- ggarrange(p1, p2, common.legend = T, legend = "bottom", nrow = 2)

pdf(paste0(main, "num_samples_and_cells_integrated.pdf"), width = 10)
print(num_cells_int)
dev.off()

rds.combined@meta.data$id_sex_age_trimester <- paste(rds.combined@meta.data$id_sex_age,rds.combined@meta.data$trimesters,  sep="_")

num_cells <- as.data.frame(table(rds.combined$id_sex_age_trimester))
num_cells <- separate(num_cells, Var1, into=c("id", "sex", "age", "trimester"), sep="_")
num_cells$age <- as.numeric(num_cells$age) 
num_cells$trimester <- factor(num_cells$trimester, c("1st", "2nd", "3rd"))

p1 <- ggplot() +
  geom_bar(data = num_cells, aes(trimester, fill=sex), position = "dodge") +
  geom_hline(yintercept = c(3), linetype = "dashed") +
  #geom_rect(data = trimesters, aes(xmin = bg_start,
  #          xmax = bg_end,
  #          ymin = - Inf,
  #          ymax = Inf,
  #          fill=colors),
  #          alpha = 0.1) +
  scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="Trimesters", y="Number of samples", fill="Sex") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

p2 <- ggplot(num_cells, aes(trimester, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Trimesters", y="Number of cells/sample", fill="Sex") +
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

num_cells_trim <- ggarrange(p1, p2, common.legend = T, legend = "bottom", nrow = 2)

pdf(paste0(main, "num_samples_and_cells_integrated_trimester.pdf"), width = 10)
print(num_cells_trim)
dev.off()

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(rds.combined) <- "RNA"

#VlnPlot(rds.combined, features = "XIST", group.by = "id_sex_age_trimester") + NoLegend()

saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")

####### CLUSTERS

DimPlot(rds.combined, reduction = "umap", group.by = "Cell.Type") # Eze
DimPlot(rds.combined, reduction = "umap", group.by = "WGCNAcluster") # Nowakowski

### Retrieve number of Eze clusters in each Cell.Type

rds.combined@meta.data[which(is.na(rds.combined@meta.data$Cell.Type)), "Cell.Type"] <- "Nowakowski"

ct_clusters <- list()
ct_names <- vector()
for (ct in unique(rds.combined@meta.data$Cell.Type)) {
  if (ct != "Nowakowski") {
    ct_names <- c(ct_names, ct)
    ct_clusters <- append(ct_clusters, list(as.numeric(unique(rds.combined@meta.data[which(rds.combined@meta.data$Cell.Type==ct), "integrated_snn_res.0.5"]))))
  }
}
names(ct_clusters) <- ct_names

Reduce(intersect, ct_clusters) 

p1 <- DimPlot(rds.combined, reduction = "umap", group.by = "integrated_snn_res.0.5", split.by = "proj")
p2 <- DimPlot(rds.combined, reduction = "umap", group.by = "Cell.Type", split.by = "proj") 
p3 <- DimPlot(rds.combined, reduction = "umap", group.by = "Cell.Type", split.by = "Cell.Type") 

eze_clusters <- ggarrange(p1, p2, nrow = 2)
pdf(paste0(main, "Eze_clusters_integrated.pdf"), width = 8)
print(eze_clusters)
dev.off()

pdf(paste0(main, "Eze_clusters_integrated_split_ct.pdf"), width = 25)
print(p3)
dev.off()


eze <- subset(rds.combined, proj=="Eze_2021")

DimPlot(eze, reduction = "umap", group.by = "Cell.Type")

Idents(eze) <- "Cell.Type"
eze_markers <- FindAllMarkers(eze, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(eze_markers, paste0(main, "Eze_markers_for_annotation.csv"))

# let's try to recluster nowa
nowa <- subset(rds.combined, proj=="Nowakowski_2017")
DefaultAssay(nowa) <- "integrated"
nowa <- FindNeighbors(nowa, dims = 1:7)
nowa <- FindClusters(nowa, resolution = 0.5)
DimPlot(nowa, reduction = "umap")

for (ct in unique(eze_markers$cluster)) {
  ct_markers <- subset(eze_markers, cluster==ct)
  top20 <- ct_markers[1:20, "gene"]
  print(DotPlot(nowa, features = unique(top20), assay = "RNA") + RotatedAxis())
}

top5 <- vector()
for (ct in unique(eze_markers$cluster)) {
  ct_markers <- subset(eze_markers, cluster==ct)
  top5 <- c(top5, ct_markers[1:5, "gene"])
}

DotPlot(nowa, features = unique(top5), assay = "RNA") + RotatedAxis()

# From this, it seems that:
# Neuronal: 5
# Mesenchymal: 3, 6
# Radial Glial: 6, 7
# Other: 6
# Neuroepithelial: 7
# IPC: maybe 1

# still not clear..

DimPlot(nowa, reduction = "umap", group.by = "WGCNAcluster")

# empty string is probably U2
nowa@meta.data[which(nowa@meta.data$WGCNAcluster== ""), "WGCNAcluster"] <- "U2"

# try to match each nowa cluster to an eze one
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

rds.combined@meta.data$cluster_final <- coalesce(rds.combined@meta.data$Cell.Type, rds.combined@meta.data$WGCNAcluster)

for (nowa_ct in names(cluster_match)) {
  rds.combined@meta.data[which(rds.combined@meta.data$WGCNAcluster==nowa_ct), "cluster_final"] <- cluster_match[[nowa_ct]]
}

rds.combined@meta.data[which(rds.combined@meta.data$cluster_final=="Nowakowski"), "cluster_final"] <- "Unknown"

DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final")

saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")


########################

rds.combined <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")
DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final")

other_cl <- subset(rds.combined, cluster_final == "Other")

DefaultAssay(other_cl) <- "integrated"
DimHeatmap(other_cl, dims = 1:15, cells = 500, balanced = TRUE)
other_cl <- JackStraw(other_cl, num.replicate = 100)
other_cl <- ScoreJackStraw(other_cl, dims = 1:20)
JackStrawPlot(other_cl, dims = 1:15)
ElbowPlot(other_cl)
other_cl <- FindNeighbors(other_cl, dims = 1:8)
other_cl <- FindClusters(other_cl, resolution = 0.5)
DimPlot(other_cl, reduction = "umap", group.by = "WGCNAcluster")


extra_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/"

cm1 <- read.csv(paste0(extra_path, "CellMarker.csv"))
cm2 <- read.csv(paste0(extra_path, "CellMarker-2.csv"))
cm3 <- read.csv(paste0(extra_path, "CellMarker-3.csv"))
cm_df <- rbind(cm1, cm2, cm3)
rm(cm1, cm2, cm3)
cm_df <- cm_df[, -c(6:8)]
col_factors <- c("Species",
                 "Tissue",
                 "Cell.Type",
                 "Cancer")
cm_df[col_factors] <- lapply(cm_df[col_factors], as.factor) 
cm_df <- subset(cm_df, subset = (Cancer == "Normal"))
cm_df <- subset(cm_df, subset = (Tissue == c("Brain") | Tissue == c("Dorsolateral prefrontal cortex")))
cm_df$Cancer <- droplevels(cm_df$Cancer)
cm_df$Tissue <- droplevels(cm_df$Tissue)
cm_df$Cell.Type <- droplevels(cm_df$Cell.Type)

ct_list <- c(
  "Astrocyte" ="astrocyte",
  "B cell" = "B",                         
  "Endothelial cell" = "EC",             
  "Glial cell" = "Glia" ,                    
  "Glutamatergic neuron" = "EN",         
  "Interstitial cell" = "IC",              
  "Lake et al.Science.Ex1" = "EN",              
  "Lake et al.Science.Ex2" = "EN",             
  "Lake et al.Science.Ex3" = "EN",              
  "Lake et al.Science.Ex4" = "EN",             
  "Lake et al.Science.Ex5" = "EN",              
  "Lake et al.Science.Ex6" = "EN",             
  "Lake et al.Science.Ex7" = "EN",              
  "Lake et al.Science.Ex8" = "EN",            
  "Lake et al.Science.In1" = "IN",              
  "Lake et al.Science.In2" = "IN",        
  "Lake et al.Science.In3" = "IN",         
  "Lake et al.Science.In4" = "IN",        
  "Lake et al.Science.In5" = "IN",         
  "Lake et al.Science.In6" = "IN",        
  "Lake et al.Science.In7" = "IN",         
  "Lake et al.Science.In8" = "IN",        
  "M1 macrophage" = "Macrophage",                  
  "M2 macrophage" = "Macrophage",                 
  "Macrophage"= "Macrophage",                      
  "Microglial cell" = "Microglia",                
  "Neural progenitor cell" = "NPC",          
  "Neural stem cell" = "NSC",                
  "Neuron" = "Neuron",                          
  "Neutrophil" = "Neutrophil",                    
  "Oligodendrocyte" = "Oligodendrocyte",                
  "Oligodendrocyte precursor cell" = "Oligodendrocyte",
  "Oligodendrocyte progenitor cell" = "Oligodendrocyte",
  "Pericyte"  = "Pericyte",                   
  "Purkinje cell" = "Neuron",                  
  "Stem cell" =    "Stem cell",                 
  "T cell"  = "T",                         
  "T helper2 (Th2) cell" = "T"
)

cm_df$ct <- rep(NA, nrow(cm_df))

for (ct in levels(cm_df$Cell.Type)) {
  cm_df[which(cm_df$Cell.Type==ct), "ct"] <- ct_list[[ct]]
}
cm_df$ct <- as.factor(cm_df$ct)
markers <- list()
for (i in levels(cm_df$ct)) {
  markers <- append(markers, list((cm_df[which(cm_df$ct==i),"Cell.Marker"])))
}
names(markers) <- levels(cm_df$ct)
for (ct in names(markers)) {
  gene_list <- vector()
  for (i in 1:length(markers[[ct]])) {
    gene_list <- c(gene_list, str_split(markers[[ct]][i], ", "))
  }
  markers[[ct]] <- unique(unlist(gene_list))
}

markers <- markers[c("astrocyte",       "B",               "EC",                       "Glia",            "IC",                     "Macrophage",     
                     "Microglia",              "Neutrophil",        "Oligodendrocyte", "Pericyte",        "Stem cell",      
                    "T" )]

for (ct in names(markers)) {
  DoHeatmap(other_cl,markers[[ct]])
}

DoHeatmap(other_cl,markers[[1]]) # -> 1,3? "astrocyte"
DoHeatmap(other_cl,markers[[2]]) # no features  "B"      
DoHeatmap(other_cl,markers[[3]]) # -> 8? "EC"
DoHeatmap(other_cl,markers[[4]]) # -> 7? "Glia"  
DoHeatmap(other_cl,markers[[5]]) # no features "IC" 
DoHeatmap(other_cl,markers[[6]]) # -> 1? "Macrophage" 
DoHeatmap(other_cl,markers[[7]]) # -> 8 "Microglia" 
DoHeatmap(other_cl,markers[[8]]) # no features "Neutrophil" 
DoHeatmap(other_cl,markers[[9]]) # -> all? "Oligodendrocyte"
DoHeatmap(other_cl,markers[[10]]) # no features "Pericyte" 
DoHeatmap(other_cl,markers[[11]]) # -> none?  "Stem cell" 
DoHeatmap(other_cl,markers[[12]]) # no features  "T"  

# it seems that only 8 can be annotated as microglia 

other_cl@meta.data$other_ann <- rep("no_data", nrow(other_cl@meta.data))
other_cl@meta.data[which(other_cl@meta.data$seurat_clusters==8), "other_ann"] <- "Microglia"

nowa_markers <- read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Nowakowski_2017_suppl.xlsx",
                         sheet = 5)

nowa_sub <- unique(other_cl$WGCNAcluster)
nowa_sub <- nowa_sub[2:length(nowa_sub)]

nowa_markers <- subset(nowa_markers, cluster %in% nowa_sub)

DoHeatmap(other_cl, unique(nowa_markers$gene)) + RotatedAxis()


saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")

#####################

rds.combined <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")
DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final")

other_cl <- subset(rds.combined, cluster_final == "Other")

DefaultAssay(other_cl) <- "integrated"
DimHeatmap(other_cl, dims = 1:15, cells = 500, balanced = TRUE)
other_cl <- JackStraw(other_cl, num.replicate = 100)
other_cl <- ScoreJackStraw(other_cl, dims = 1:20)
JackStrawPlot(other_cl, dims = 1:15)
ElbowPlot(other_cl)
other_cl <- FindNeighbors(other_cl, dims = 1:8)
other_cl <- FindClusters(other_cl, resolution = 0.5)
DimPlot(other_cl, reduction = "umap", group.by = "WGCNAcluster")


expr_mat_all <- GetAssayData(other_cl[["RNA"]], slot="data")

other_cl@meta.data[which(is.na(other_cl@meta.data$WGCNAcluster)), "WGCNAcluster"] <- "Eze"

other_cl@meta.data$seurat_cluster_combo <- paste(other_cl@meta.data$seurat_clusters, other_cl@meta.data$WGCNAcluster, sep = "_")

Idents(other_cl) <- "seurat_cluster_combo"

cell_info <- data.frame()
for (i in unique(other_cl@meta.data$seurat_cluster_combo)) {
  print(i)
  cell_id <- WhichCells(other_cl, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}
cell_info <- separate(cell_info, og_group, into = c("seurat_cluster", "WGCNAcluster"), remove = F, sep = "_")

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
 
#saveRDS(expr_mat_all, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/top_2000_SD_expr_matrix.rds")


avg_expr <- data.frame(expr_mat_all[,1])
col_names <- vector()

for (cl_seu in unique(cell_info$seurat_cluster)) {
  seu_cells <- cell_info[which(cell_info$seurat_cluster==cl_seu), "cell_id"]
  avg_expr <- cbind(avg_expr, rowMeans(expr_mat_all[,seu_cells]))
  col_names <- c(col_names, paste0("seurat_", cl_seu))
}
for (cl_WGCNA in unique(cell_info$WGCNAcluster)) {
  WGCNA_cells <- cell_info[which(cell_info$WGCNAcluster==cl_WGCNA), "cell_id"]
  avg_expr <- cbind(avg_expr, rowMeans(expr_mat_all[,WGCNA_cells]))
  col_names <- c(col_names, cl_WGCNA)
}
colnames(avg_expr) <- c("Genes", col_names)
rownames(avg_expr) <- avg_expr$Genes
avg_expr[,1] <- NULL


cor_pearson <- cor(avg_expr, method = "pearson")
p.mat <- cor_pmat(avg_expr)

library(ggcorrplot)
ggcorrplot(cor_spearman, p.mat = p.mat, type = "lower", lab = T)
ggsave(paste0(main, "Eze_Nowakowski_integrated/other_cluster_correlation_mtx.pdf"))

# the only somewhat stronger correlation is between Seurat_8 and microglia, which is what we had observed previously

other_cl@meta.data$other_ann <- rep(NA, nrow(other_cl@meta.data))
other_cl@meta.data[which(other_cl@meta.data$seurat_clusters==8), "other_ann"] <- "Microglia"

other_cl@meta.data$other_ann <- coalesce(other_cl@meta.data$other_ann, other_cl@meta.data$WGCNAcluster)

DimPlot(other_cl, reduction = "umap",  group.by = "other_ann")

Idents(other_cl) <- "other_ann"
microglia_cells <- WhichCells(other_cl, idents = "Microglia")

rds.combined@meta.data[which(rds.combined@meta.data$Cell %in% microglia_cells), "cluster_final"] <- "Microglia"


saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")

rds.combined@meta.data$id_sex_age_cluster <- paste(rds.combined@meta.data$id_sex_age, rds.combined@meta.data$cluster_final, sep="_")
rds.combined@meta.data$proj_id_sex_age_cluster <- paste(rds.combined@meta.data$proj, rds.combined@meta.data$id_sex_age_cluster, sep="_")

rds.combined@meta.data$id_sex_trim_cluster <- paste(rds.combined@meta.data$id_sex_age_trimester, rds.combined@meta.data$cluster_final, sep="_")
rds.combined@meta.data$proj_id_sex_trim_cluster <- paste(rds.combined@meta.data$proj, rds.combined@meta.data$id_sex_trim_cluster, sep="_")


num_cells <- as.data.frame(table(rds.combined$proj_id_sex_trim_cluster))
num_cells <- separate(num_cells, Var1, into=c("proj","year", "id", "sex", "age", "trim", "cluster"), sep="_")
num_cells$proj <- paste(num_cells$proj, num_cells$year, sep="_")
num_cells$year <- NULL
num_cells$age <- as.numeric(num_cells$age) 

num_cells_sex_ct <- ggplot(num_cells, aes(cluster, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  labs(x="Cell types", y="Number of cells/sample", fill="Sex") +
  facet_wrap(~trim, scales = "free") +
  geom_hline(yintercept = 500, linetype="dashed") +
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

pdf(paste0(main, "Eze_Nowakowski_integrated/num_cells_ct_trim_sex.pdf"))
print(num_cells_sex_ct)
dev.off()




###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################


# Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html
input_rds_path <-  "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC"
input_rds_files <- list.files(path = input_rds_path, pattern = ".rds", full.names = T)[c(3,8)]
input_rds <- lapply(input_rds_files,function(x) {
  readRDS(file = x)
})
names(input_rds) <- list.files(path = input_rds_path, pattern = ".rds", full.names = F)[c(3,8)]
names(input_rds) <- str_remove_all(names(input_rds), ".rds")

# normalize and identify variable features for each dataset independently
input_rds <- lapply(X = input_rds, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = input_rds)
common.anchors <- FindIntegrationAnchors(object.list = input_rds, anchor.features = features)

# this command creates an 'integrated' data assay
rds2.combined <- IntegrateData(anchorset = common.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(rds2.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
rds2.combined <- ScaleData(rds2.combined, verbose = FALSE)
rds2.combined <- RunPCA(rds2.combined, npcs = 30, verbose = FALSE)

DimHeatmap(rds2.combined, dims = 1:15, cells = 500, balanced = TRUE)

rds2.combined <- RunUMAP(rds2.combined, reduction = "pca", dims = 1:15)
rds2.combined <- FindNeighbors(rds2.combined, reduction = "pca", dims = 1:15)
rds2.combined <- FindClusters(rds2.combined, resolution = 0.5)

rds2.combined@project.name <- "Nowakowski_Velmeshev_integrated"

# Visualization
DimPlot(rds2.combined, reduction = "umap")
DimPlot(rds2.combined, reduction = "umap", group.by = "proj")
ggsave(paste0(main, rds2.combined@project.name, "_proj.pdf"))

rds2.combined@meta.data$age_final <- coalesce(as.character(rds2.combined@meta.data$Age_in_Weeks), rds2.combined@meta.data$age)

rds2.combined@meta.data$trimester <- rep("3rd", nrow(rds2.combined@meta.data))

rds2.combined@meta.data[which(0 < rds2.combined@meta.data$Age_in_Weeks & rds2.combined@meta.data$Age_in_Weeks < 13), "trimester"] <- "1st"
rds2.combined@meta.data[which(12 < rds2.combined@meta.data$Age_in_Weeks & rds2.combined@meta.data$Age_in_Weeks < 27), "trimester"] <- "2nd"
rds2.combined@meta.data[which(26 < rds2.combined@meta.data$Age_in_Weeks & rds2.combined@meta.data$Age_in_Weeks < 41), "trimester"] <- "3rd"

DimPlot(rds2.combined, reduction = "umap", group.by = "trimester", split.by = "proj")
ggsave(paste0(main, rds2.combined@project.name, "_trimester_per_proj.pdf"))

DimPlot(rds2.combined, reduction = "umap", group.by = "WGCNAcluster", split.by = "trimester")
ggsave(paste0(main, rds2.combined@project.name, "_trimester_WGCNAcluster.pdf"))


saveRDS(rds2.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Nowakowski_Velmeshev_3rd_trimester_integrated.rds")
#rm(input_rds, input_rds_files, input_rds_path)

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################



# Their tutorial - https://satijalab.org/seurat/articles/integration_introduction.html
# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

# normalize and identify variable features for each dataset independently
ifnb.list <- lapply(X = ifnb.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = ifnb.list)
immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, anchor.features = features)

# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(immune.combined) <- "RNA"
nk.markers <- FindConservedMarkers(immune.combined, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
head(nk.markers)

FeaturePlot(immune.combined, features = c("CD3D", "SELL", "CREM", "CD8A", "GNLY", "CD79A", "FCGR3A",
                                          "CCL2", "PPBP"), min.cutoff = "q9")


immune.combined <- RenameIdents(immune.combined, `0` = "CD14 Mono", `1` = "CD4 Naive T", `2` = "CD4 Memory T",
                                `3` = "CD16 Mono", `4` = "B", `5` = "CD8 T", `6` = "NK", `7` = "T activated", `8` = "DC", `9` = "B Activated",
                                `10` = "Mk", `11` = "pDC", `12` = "Eryth", `13` = "Mono/Mk Doublets", `14` = "HSPC")
DimPlot(immune.combined, label = TRUE)

Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                      "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                      "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()

library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
t.cells <- subset(immune.combined, idents = "CD4 Naive T")
Idents(t.cells) <- "stim"
avg.t.cells <- as.data.frame(log1p(AverageExpression(t.cells, verbose = FALSE)$RNA))
avg.t.cells$gene <- rownames(avg.t.cells)

cd14.mono <- subset(immune.combined, idents = "CD14 Mono")
Idents(cd14.mono) <- "stim"
avg.cd14.mono <- as.data.frame(log1p(AverageExpression(cd14.mono, verbose = FALSE)$RNA))
avg.cd14.mono$gene <- rownames(avg.cd14.mono)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
b.interferon.response <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", verbose = FALSE)
head(b.interferon.response, n = 15)

FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3,
            cols = c("grey", "red"))

plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

