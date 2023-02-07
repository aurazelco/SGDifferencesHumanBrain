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
`%!in%` <- Negate(`%in%`)

main_outs <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Eze_Nowakowski_integrated/"
main_deg <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Eze_Nowakowski_integrated_2nd_trimester/"
main_scenic <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/Eze_Nowakowski_integrated_2nd_trimester/"
input_rds_path <-  "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC"

####### Modified tutorial from https://satijalab.org/seurat/articles/integration_introduction.html
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
ggsave(paste0(main_outs, rds.combined@project.name, "_proj.pdf"))

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
ggsave(paste0(main_outs, rds.combined@project.name, "_age_combined.pdf"))

# 1st_trimester: 1 - 12 
# 2nd_trimester: 13 - 26
# 3rd_trimester: 27 - 40

rds.combined@meta.data$trimesters <- rep("no_data", nrow(rds.combined@meta.data))

rds.combined@meta.data[which(0 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 13), "trimesters"] <- "1st"
rds.combined@meta.data[which(12 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 27), "trimesters"] <- "2nd"
rds.combined@meta.data[which(26 < rds.combined@meta.data$age_final & rds.combined@meta.data$age_final < 41), "trimesters"] <- "3rd"

DimPlot(rds.combined, reduction = "umap", group.by = "trimesters")
ggsave(paste0(main_outs, rds.combined@project.name, "_trimesters.pdf"))


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

pdf(paste0(main_outs, "num_samples_and_cells_integrated.pdf"), width = 10)
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

pdf(paste0(main_outs, "num_samples_and_cells_integrated_trimester.pdf"), width = 10)
print(num_cells_trim)
dev.off()

# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(rds.combined) <- "RNA"

VlnPlot(rds.combined, features = "XIST", group.by = "id_sex_age_trimester") + NoLegend()

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
pdf(paste0(main_outs, "Eze_clusters_integrated.pdf"), width = 8)
print(eze_clusters)
dev.off()

pdf(paste0(main_outs, "Eze_clusters_integrated_split_ct.pdf"), width = 25)
print(p3)
dev.off()


eze <- subset(rds.combined, proj=="Eze_2021")

DimPlot(eze, reduction = "umap", group.by = "Cell.Type")

Idents(eze) <- "Cell.Type"
eze_markers <- FindAllMarkers(eze, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.1)
write.csv(eze_markers, paste0(main_outs, "Eze_markers_for_annotation.csv"))

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
DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final", split.by = "sex")
ggsave(paste0(main_outs, "clusters_by_sex.pdf"))

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
ggsave(paste0(main_outs, "other_cluster_correlation_mtx.pdf"))

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

pdf(paste0(main_outs, "num_cells_ct_trim_sex.pdf"))
print(num_cells_sex_ct)
dev.off()

########################

rds.combined <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated.rds")
DimPlot(rds.combined, reduction = "umap", group.by = "cluster_final")
DimPlot(rds.combined, reduction = "umap", group.by = "trimesters")


trim_2nd <- subset(rds.combined, trimesters=="2nd")
rm(rds.combined)


DefaultAssay(trim_2nd) <- "integrated"
DimHeatmap(trim_2nd, dims = 1:15, cells = 500, balanced = TRUE)
trim_2nd <- JackStraw(trim_2nd, num.replicate = 100)
trim_2nd <- ScoreJackStraw(trim_2nd, dims = 1:20)
JackStrawPlot(trim_2nd, dims = 1:15)
ElbowPlot(trim_2nd)
trim_2nd <- FindNeighbors(trim_2nd, dims = 1:13)
trim_2nd <- FindClusters(trim_2nd, resolution = 0.5)
DimPlot(trim_2nd, reduction = "umap", group.by = "cluster_final")
DimPlot(trim_2nd, reduction = "umap", group.by = "trimesters")
saveRDS(trim_2nd, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated_2nd_trimester.rds")

########################

trim_2nd <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/Eze_Nowakowski_integrated_2nd_trimester.rds")
trim_2nd@meta.data$sex_ct <- paste(trim_2nd@meta.data$sex, trim_2nd@meta.data$cluster_final, sep="_")

trim_2nd@project.name <- "Eze_Nowakowski_integrated_2nd_trimester"

num_sex_ct <- as.data.frame(table(trim_2nd$sex_ct))
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

trim_2nd_output <- paste0(main_deg, trim_2nd@project.name, "/outputs")

dir.create(trim_2nd_output, recursive = T, showWarnings = F)

for (min_cells in min_num_cells) {
  num_filt <- FiltDF(num_sex_ct, min_cells)
  write.csv(num_filt, file = paste0(trim_2nd_output, "/final_filt_", min_cells, ".csv"),
            row.names = F)
  pdf(paste0(trim_2nd_output, "/filt_counts_", min_cells, ".pdf"), 10, 15)
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

saveRDS(trim_2nd, paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))


############ For 02C_Conservation

trim_2nd <- readRDS(paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

Idents(trim_2nd) <- "sex_ct"

expr_mat_all_cts <- GetAssayData(trim_2nd[["RNA"]], slot="data")

cell_info <- data.frame()
for (i in unique(trim_2nd@meta.data$sex_ct)) {
  print(i)
  cell_id <- WhichCells(trim_2nd, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}

write.csv(cell_info, paste0(main_deg, "cell_info.csv"))

rm(trim_2nd)
expr_mat_all_cts <- as.data.frame(as.matrix(expr_mat_all_cts))

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info$og_group)) {
  og_cells <- c("Genes", cell_info[which(cell_info$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all_cts[ , (names(expr_mat_all_cts) %in% og_cells)]
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

write.csv(tot_genes, paste0(main_deg, "tot_genes_ct.csv"))

############ For SCENIC - HIGHEST VARIABLE GENES

trim_2nd <- readRDS(paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

trim_2nd@meta.data$sample_id <- rownames(trim_2nd@meta.data)

trim_2nd@meta.data$sex_ct_sample <- paste(trim_2nd@meta.data$sex_ct, trim_2nd@meta.data$sample_id, sep="_")

Idents(trim_2nd) <- "sex_ct_sample"

expr_mat_all <- GetAssayData(trim_2nd[["RNA"]], slot="data")

cell_info <- read.csv(paste0(main_deg, "cell_info.csv"))
cell_info$X <- NULL

saveRDS(trim_2nd, paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

rm(trim_2nd)

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

dir.create(main_scenic, showWarnings = F, recursive = T)

saveRDS(expr_mat_all, paste0(main_scenic, "top_2000_SD_expr_matrix.rds"))

##### Map the samples back to the groups they belong to

expr_mat_all <- readRDS(paste0(main_scenic, "top_2000_SD_expr_matrix.rds"))

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

group_list100 <- remove_dfs(group_list, 100)
group_list500 <- remove_dfs(group_list, 500)

plot_group_numbers <- function(df_list, thresh) {
  ids <- as.data.frame(names(df_list))
  colnames(ids) <- c("groups")
  ids <- separate(ids, groups, into = c("sex","ct"), sep="_", remove=FALSE)
  col_factors <- c("sex","ct")
  ids[col_factors] <- lapply(ids[col_factors], as.factor)  
  ids$length_groups <- sapply(1:length(names(df_list)), function(i) ncol(df_list[[i]]))
  pdf(paste0(main_scenic, "num_filt_", thresh, "_cells.pdf"))
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

expr_mat_all <- readRDS(paste0(main_scenic, "top_2000_SD_expr_matrix.rds"))
cell_info <- read.csv(paste0(main_deg, "cell_info.csv"))
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
  names(sampled_dfs) <- lapply(1:length(sampled_names), function(i) str_replace_all(sampled_names[i]," ", "-"))
  dir.create(paste0(main, "sampled_", num_cells, "_cells"), showWarnings = FALSE)
  lapply(1:length(names(sampled_dfs)), function(i) write.csv(sampled_dfs[[i]], 
                                                            file = paste0(main, "sampled_", num_cells, "_cells/", names(sampled_dfs)[i], ".csv"),
                                                            row.names = FALSE))
  return(sampled_dfs)
}


df_100_sampled <- rand_sample(df_list100, 3, 100, main_scenic)

############ For SCENIC - HIGHEST VARIABLE GENES

trim_2nd <- readRDS(paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

trim_2nd@meta.data$id_sex <- paste(trim_2nd@meta.data$id, trim_2nd@meta.data$sex, sep="_")

VlnPlot(trim_2nd, features = "XIST", group.by = "id_sex", assay = "RNA") + NoLegend()
ggsave(paste0(main_outs, "XIST_RNA_expression.pdf"))

############ NUM CELLS

trim_2nd <- readRDS(paste0(input_rds_path, "/Eze_Nowakowski_integrated_2nd_trimester.rds"))

trim_2nd@meta.data$proj_sex_ct <- paste(trim_2nd@meta.data$proj, trim_2nd@meta.data$sex,trim_2nd@meta.data$cluster_final, sep="/")

num_cells <- as.data.frame(table(trim_2nd@meta.data$proj_sex_ct))
num_cells <- separate(num_cells, Var1, into = c("project", "sex", "ct"), sep = "/")
write.csv(num_cells, paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Eze_Nowakowski_num_cells.csv"))