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

DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final")

saveRDS(rds.combined, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")


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

