setwd("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/")
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)

# load the dataset
disco_brain <- readRDS("brainV1.0.rds")

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

# cause_of_death, age_group, age, bmi throw error because of NAs

length(unique(disco_brain$sample))

allmisscols <- sapply(disco_brain@meta.data, function(x) all(is.na(x) | x == '' ))
allmisscols # -> these columns have NAs, so cannot be plotted

#disco_brain_noNA <- subset(disco_brain, subset = age != "NA") 
# this removes NAs, but the age groups are not including neonates -> 
# GSE134355 and GSE165388 from manually parsed dataset

#disco_fetal <- subset(disco_brain, subset = project_id == c("GSE134355","GSE165388")) 
# these projects are not included in the RDS file for some reason...

# extracted the sex NAs to see if we can assign it based on XIST expression and X genes from O'Brien

FeaturePlot(disco_brain, features = "XIST", split.by = "gender")

sex_NA <- subset(disco_brain, subset = gender !=  c("F","M") ) # trying to subset based on NA did not work


DimPlot(sex_NA, reduction = "umap", group.by = "gender") # NA


library(readxl)
Xgenes <- read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/O'Brien_bioRXiv_suppl.xlsx",
                     sheet="TblS4",
                     col_names = TRUE,
                     skip=2)

length(unique(sex_NA$sample)) # 103

FeaturePlot(sex_NA, features = "XIST")
VlnPlot(sex_NA, features = "XIST", group.by = "project_id")
VlnPlot(sex_NA, features = "XIST", group.by = "sample", ncol=10)


xist <- AverageExpression(sex_NA, features = "XIST", group.by = "project_id")
xist_samples <- AverageExpression(sex_NA, features = "XIST", group.by = "sample")

xist_df <- as.data.frame(xist[2])
rows <- colnames(xist_df)
xist_df <- transpose(xist_df)
rownames(xist_df) <- rows
colnames(xist_df) <- c("XIST")
setDT(xist_df, keep.rownames = "proj")

xist_sample_df <- as.data.frame(xist_samples[2])
rows <- colnames(xist_sample_df)
xist_sample_df <- transpose(xist_sample_df)
rownames(xist_sample_df) <- rows
colnames(xist_sample_df) <- c("XIST")
setDT(xist_sample_df, keep.rownames = "samples")

ggplot(xist_df, aes(proj,XIST,fill=proj)) +
  geom_bar(stat="identity") +   
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))

ggplot(xist_sample_df, aes(samples,XIST,fill=samples)) +
  geom_bar(stat="identity", show.legend = FALSE) +   
  geom_hline(yintercept=1, linetype="dashed", color = "black") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"))

# Violin and FeaturePlot agree with expression levels, but I do not get the same numbers if I extract the average expression

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5816366/
# not found: "AZFa","TB4Y",  "AZFb","AZFc", "XKRY", "VCY1",,  "VCY2"
# not expressed: "TGIF2LY", "TSPY1", "CYorf15", "RPS4Y2","PRY", "PRY2", "RBMY1A1", "DAZ1", "DAZ2", "DAZ3", "DAZ4", "CDY1,
# CDY2",
Ygenes <- c("SRY", "ZFY", "RPS4Y1", "AMELY", "TBL1Y", "PCDH11Y", "TSPY2",  
            "USP9Y", "DDX3Y", "UTY",  "EIF1AY", "KDM5D", 
            "HSFY1", "HSFY2")

for (gene in Ygenes) {
  FeaturePlot(sex_NA, features = gene)
}

Ygenes_filt <- c("ZFY", "RPS4Y1", "PCDH11Y",  
                 "USP9Y", "DDX3Y", "UTY",  
                 "EIF1AY", "KDM5D")

FeaturePlot(sex_NA, features = Ygenes)
FeaturePlot(sex_NA, features = Ygenes_filt)
FeaturePlot(disco_brain, features = Ygenes_filt)

# chose DDX3Y as "male" marker

# -> subset again, using the FeaturePlot limit - totally need to check if 2.5 is OK
female_NA <- subset(sex_NA, subset = XIST >= 2.5)
male_NA <- subset(sex_NA, subset = DDX3Y >= 2.5)

females <- unique(female_NA$sample)
males <- unique(male_NA$sample)

length(intersect(females, males)) # 17 -> samples found in both
length(setdiff(females, males)) # 27 -> only in females
length(setdiff(males, females)) # 59 -> only in males

mix <- intersect(females, males)

mix_NA <- subset(sex_NA, subset = sample == c(mix))


FeaturePlot(female_NA, features = "DDX3Y")
FeaturePlot(male_NA, features = "XIST")

# not sure how to proceed... maybe I should take in consideration only DDX3Y and 
# label everything else F?

# maybe for the moment I will focus on the known ones

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

ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/analysis/DISCOv1.0/DISCO_FM/samples_counts.png")

FeaturePlot(disco_FM, features = "XIST", split.by = "gender", reduction = "umap")

ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/analysis/DISCOv1.0/DISCO_FM/xist_exp.png")

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

write.csv(celltype_num, "celltype_num.csv")

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

ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/analysis/DISCOv1.0/DISCO_FM/ct_counts.png")

saveRDS(disco_FM, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DISCOv1.0_FM/disco_FM.rds")

########### August 11th

disco <- readRDS("brainV1.0_all_FM.rds")

celltype_num <- read.csv("celltype_num.csv")
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

write.csv(celltype_num, "celltype_num_filt.csv")

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

saveRDS(disco_filt, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

Idents(disco_filt) <- "proj_sex_ct"

DimPlot(disco_filt, reduction = "umap")

GSE157827_sex_microglia <- FindMarkers(disco_filt, 
                                      ident.1 = "GSE157827_F_microglia", 
                                      ident.2 = "GSE157827_M_microglia",
                                      logfc.threshold = 0.25,
                                      min.pct = 0.1)

GSE153807_sex_microglia <- FindMarkers(disco_filt, 
                                       ident.1 = "GSE153807_F_microglia", 
                                       ident.2 = "GSE153807_M_microglia",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1)

PRJNA544731_sex_microglia <- FindMarkers(disco_filt, 
                                       ident.1 = "PRJNA544731_F_microglia", 
                                       ident.2 = "PRJNA544731_M_microglia",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1)

GSE174367_sex_microglia <- FindMarkers(disco_filt, 
                                       ident.1 = "GSE174367_F_microglia", 
                                       ident.2 = "GSE174367_M_microglia",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1)


write.csv(GSE157827_sex_microglia, "GSE157827_sex_microglia.csv")
write.csv(GSE153807_sex_microglia, "GSE153807_sex_microglia.csv")
write.csv(PRJNA544731_sex_microglia, "PRJNA544731_sex_microglia.csv")
write.csv(GSE174367_sex_microglia, "GSE174367_sex_microglia.csv")


disco_filt@meta.data$sex_ct <- paste(disco_filt@meta.data$gender, disco_filt@meta.data$ct, sep="_")

Idents(disco_filt) <- "sex_ct"

sex_microglia <- FindMarkers(disco_filt, 
                            ident.1 = "F_microglia", 
                            ident.2 = "M_microglia",
                            logfc.threshold = 0.25,
                            min.pct = 0.1)


write.csv(sex_microglia, "overall_sex_microglia.csv")


sex_microglia <- read.csv("overall_sex_microglia.csv", row.names = 1)

top10_overall <- rownames(sex_microglia)[1:10]

VlnPlot(disco_microglia, features = top10_overall, group.by = "gender")
RidgePlot(disco_microglia, features = top10_overall, group.by = "gender")

# Number of cells per analysis

overall_microglia_num <- as.data.frame(table(disco_filt$sex_ct))

proj_microglia_num <- as.data.frame(table(disco_filt$proj_sex_ct))

numbers_df <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses.csv")
numbers_df[,1] <- NULL

label <- c("DISCO_microglia")
f <- overall_microglia_num[which(overall_microglia_num$Var1=="F_microglia"), "Freq"]
m <- overall_microglia_num[which(overall_microglia_num$Var1=="M_microglia"), "Freq"]

proj <- unique(disco_filt$project_id)

for (id in proj) {
  id_label <- paste(id, "microglia", sep="_")
  label <- c(label, id_label)
  groupF <- paste(id, "F_microglia", sep="_")
  f <- c(f, proj_microglia_num[which(proj_microglia_num$Var1==groupF), "Freq"])
  groupM <- paste(id, "M_microglia", sep="_")
  m <- c(m, proj_microglia_num[which(proj_microglia_num$Var1==groupM), "Freq"])
}

numbers_df <- rbind(numbers_df, data.frame(label, f, m))
colnames(numbers_df) <- c("Analysis", "F", "M")
numbers_df$Analysis <- as.factor(numbers_df$Analysis)
numbers_df$dataset <- rep("DISCO", length.out=length(numbers_df$Analysis))
numbers_df[c(1:2), "dataset"] <- "Allen"
numbers_df$dataset <- as.factor(numbers_df$dataset)

write.csv(numbers_df, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses_all.csv")


library(data.table) # to re-arrange dataframes

df1 <- melt(as.data.table(numbers_df, "Analysis"))
df1[,1] <- NULL
colnames(df1) <- c("analysis", "dataset", "sex", "count")
df1$dataset <- as.factor(df1$dataset)
df1$sex <- as.factor(df1$sex)



ggplot(df1, aes(analysis, count, fill=sex)) +
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="Groups", y="Nuclei count", fill="Sex") +
  facet_wrap(~dataset, scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")
  
ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses.jpeg", 
       dpi=300)


disease <- as.data.frame(table(disco_filt$disease))


disco_filt@meta.data$proj_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data$proj_sex_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, sep="_")

DimPlot(disco_filt, reduction = "umap", group.by = "proj_sex_disease")


disease <- as.data.frame(table(disco_filt$proj_sex_disease))
disease <- separate(disease, Var1, into = c("project", "sex", "disease"), sep = "_")

col_factors <- c("project", 
                 "sex",
                 "disease"
)
disease[col_factors] <- lapply(disease[col_factors], as.factor)  
names(disease)[names(disease) == 'Freq'] <- "count"

write.csv(disease, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses_disease.csv")


df2 <- melt(as.data.table(disease, "project"))
df2[,1] <- NULL
df2[,4] <- NULL
colnames(df2) <- c("project", "sex", "disease", "count")

ggplot(df2, aes(disease, count, fill=sex)) +
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="Disease status", y="Nuclei count", fill="Sex") +
  facet_wrap(~project, scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses_disease.jpeg", 
       dpi=300)

FeaturePlot(disco_filt, reduction = "umap", features= "TYROBP")
DimPlot(disco_filt, reduction = "umap", group.by = "ct")

FeaturePlot(disco_microglia, reduction = "umap", features= "TYROBP")
DimPlot(disco_mi, reduction = "umap", group.by = "ct")

saveRDS(disco_filt, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

disco_filt@meta.data[which(is.na(disco_filt@meta.data$disease)), "disease"] <- "Unknown"

unique(disco_filt@meta.data[which(disco_filt@meta.data$disease=="Unknown"), "sample_type"])

disco_filt@meta.data[which(disco_filt@meta.data$disease=="Unknown"), "disease"] <- "Normal"

disco_filt@meta.data$proj_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$disease, sep="_")
disco_filt@meta.data$proj_sex_disease <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, sep="_")

disease <- as.data.frame(table(disco_filt$proj_sex_disease))
disease <- separate(disease, Var1, into = c("project", "sex", "disease"), sep = "_")

col_factors <- c("project", 
                 "sex",
                 "disease"
)
disease[col_factors] <- lapply(disease[col_factors], as.factor)  
names(disease)[names(disease) == 'Freq'] <- "count"

write.csv(disease, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses_disease.csv")

df2 <- melt(as.data.table(disease, "project"))
df2[,1] <- NULL
df2[,4] <- NULL
colnames(df2) <- c("project", "sex", "disease", "count")

ggplot(df2, aes(disease, count, fill=sex)) +
  geom_bar(stat="identity", position = "dodge") + 
  labs(x="Disease status", y="Nuclei count", fill="Sex") +
  facet_wrap(~project, scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom")

ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/number_analyses_disease.jpeg", 
       dpi=300)

disco_filt@meta.data$proj_sex_disease_ct <- paste(disco_filt@meta.data$proj_sex_disease, disco_filt@meta.data$ct, sep="_")




# DEG analysis -> only normal samples

Idents(disco_filt) <- "proj_sex_disease_ct"

disco_microglia <- subset(disco_filt, subset = ct == "microglia")

disco_microglia <- subset(disco_microglia, subset = sample_type == "normal")

normal_microglia <- as.data.frame(table(disco_microglia$proj_sex_disease_ct))

# only two projects have both M and F normal microglia

GSE157827_normal_F_microglia <- FindMarkers(disco_filt, 
                                       ident.1 = "GSE157827_F_Normal_microglia", 
                                       ident.2 = "GSE157827_M_Normal_microglia",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1,
                                       only.pos = TRUE)

GSE174367_normal_F_microglia <- FindMarkers(disco_filt, 
                                       ident.1 = "GSE174367_F_Normal_microglia", 
                                       ident.2 = "GSE174367_M_Normal_microglia",
                                       logfc.threshold = 0.25,
                                       min.pct = 0.1,
                                       only.pos = TRUE)


write.csv(GSE157827_normal_F_microglia, "20220815/GSE157827_normal_F_microglia.csv")
write.csv(GSE174367_normal_F_microglia, "20220815/GSE174367_normal_F_microglia.csv")

GSE157827_normal_M_microglia <- FindMarkers(disco_filt, 
                                            ident.1 = "GSE157827_M_Normal_microglia", 
                                            ident.2 = "GSE157827_F_Normal_microglia",
                                            logfc.threshold = 0.25,
                                            min.pct = 0.1,
                                            only.pos = TRUE)

GSE174367_normal_M_microglia <- FindMarkers(disco_filt, 
                                            ident.1 = "GSE174367_M_Normal_microglia", 
                                            ident.2 = "GSE174367_F_Normal_microglia",
                                            logfc.threshold = 0.25,
                                            min.pct = 0.1,
                                            only.pos = TRUE)


write.csv(GSE157827_normal_M_microglia, "20220815/GSE157827_normal_M_microglia.csv")
write.csv(GSE174367_normal_M_microglia, "20220815/GSE174367_normal_M_microglia.csv")


disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

num_proj_sex_disease_ct <- as.data.frame(table(disco_filt$proj_sex_disease_ct))
num_proj_sex_disease_ct <- separate(num_proj_sex_disease_ct, Var1, into = c("proj", "sex" , "disease", "ct"), sep = "_")

print(subset(num_proj_sex_disease_ct, subset = ct == "microglia" & proj %in% c("GSE157827", "GSE174367") & disease == "Normal"))

saveRDS(disco_filt, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

disco_filt@meta.data[which(disco_filt@meta.data$disease=="Alzheimer’s disease"), "disease"] <- "Alzheimer's disease"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="multiple sclerosis"), "disease"] <- "Multiple Sclerosis"  
disco_filt@meta.data[which(disco_filt@meta.data$disease=="neurosurgical disease"), "disease"] <- "Neurosurgical Disease"  

disco_filt@meta.data$proj_sex_disease_ct <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$gender, disco_filt@meta.data$disease, disco_filt@meta.data$ct, sep="_")

disco_filt <- subset(disco_filt, project_id != "GSE153807")

num_proj_sex_disease_ct <- as.data.frame(table(disco_filt$proj_sex_disease_ct))
num_proj_sex_disease_ct <- separate(num_proj_sex_disease_ct, Var1, into = c("proj", "sex" , "disease", "ct"), sep = "_")

col_factors <- c("proj", 
                 "sex",
                 "disease",
                 "ct"
)
num_proj_sex_disease_ct[col_factors] <- lapply(num_proj_sex_disease_ct[col_factors], as.factor)  
names(num_proj_sex_disease_ct)[names(num_proj_sex_disease_ct) == 'Freq'] <- "count"

for (id in levels(num_proj_sex_disease_ct$proj)) {
  pdf(paste0("20220817/outputs/", id, "_counts.pdf"), 10, 15)
  print(ggplot(num_proj_sex_disease_ct[which(num_proj_sex_disease_ct$proj==id),], aes(disease, count, fill=sex)) +
          geom_bar(stat="identity", position = "dodge") + 
          labs(x="", y="Nuclei count", fill="Sex") +
          facet_wrap(~ct, scales = "free") +
          geom_hline(yintercept = 10, linetype="dashed") +
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

FiltDF <- function(df, disease) {
  `%!in%` <- Negate(`%in%`)
  df <- droplevels(df)
  incomplete_ct <- vector()
  for (type in levels(df$ct)) {
    if ((nrow(subset(df, subset = ct == type))%%2!=0) | (any(subset(df, subset = ct == type)[,5]<10))) {
      #print(subset(df, subset = ct == type))
      incomplete_ct <- c(incomplete_ct, type)
    } else if ((nrow(subset(df, subset = ct == type))==2) & (disease!="MS")) {
      incomplete_ct <- c(incomplete_ct, type)
    }
  }
  incomplete_df <- df[df$ct %in% incomplete_ct,]
  incomplete_df <- droplevels(incomplete_df)
  incomplete_proj <- vector()
  for (type in levels(incomplete_df$ct)) {
    for (id in levels(incomplete_df$proj)) {
      if ((nrow(subset(df, subset = (ct==type & proj==id)))%%2!=0) | (any(subset(df, subset = (ct==type & proj==id))[,5]<10))) {
        incomplete_proj <- c(incomplete_proj, (paste(id, type, sep="_")))
      }
    }
  }
  df$og <- paste(df$proj, df$ct, sep="_")
  df_filt <- df[df$og %!in% incomplete_proj,]
  return(df_filt)
}

num_normal <- subset(num_proj_sex_disease_ct, disease == "Normal")
num_normal_filt <- FiltDF(num_normal, "Normal")

num_AD <- subset(num_proj_sex_disease_ct, disease == "Alzheimer's disease")
num_AD_filt <- FiltDF(num_AD, "AD")

num_MS <- subset(num_proj_sex_disease_ct, disease == "Multiple Sclerosis")
num_MS_filt <- FiltDF(num_MS, "MS")

#num_list <- list(num_normal_filt, num_AD_filt, num_MS)
#names(num_list) <- c("num_normal_filt", "num_AD_filt", "num_MS_filt")
#lapply(1:length(num_list), function(i) write.csv(num_list[i], 
#                                                 file = paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0//20220817_DEGs/outputs/",
#                                                               names(num_list[i]), ".csv"),
#                                                 row.names = TRUE))


num_filt <- rbind(num_normal_filt, num_AD_filt, num_MS_filt)
num_filt$idents <- paste(num_filt$proj, num_filt$sex, num_filt$disease, num_filt$ct, sep="_")

write.csv(num_filt, file = paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/", 
                                  "/20220817_DEGs/outputs/final_filt.csv"),
          row.names = TRUE)

for (id in levels(num_filt$proj)) {
  pdf(paste0("20220817/outputs/", id, "_filt_counts.pdf"), 10, 15)
  print(ggplot(num_filt[which(num_filt$proj==id),], aes(disease, count, fill=sex)) +
          geom_bar(stat="identity", position = "dodge") + 
          labs(x="", y="Nuclei count", fill="Sex") +
          facet_wrap(~ct, scales = "free") +
          geom_hline(yintercept = 10, linetype="dashed") +
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


saveRDS(disco_filt, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")


# extract all DEGs among different cts, for conservation plots
disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

#disco_filt@meta.data$proj_disease_ct <- paste(disco_filt@meta.data$project_id, disco_filt@meta.data$disease, disco_filt@meta.data$ct, sep="_")
Idents(disco_filt) <- "ct"

num_proj_disease_ct <- as.data.frame(table(disco_filt$proj_disease_ct))
num_proj_disease_ct <- separate(num_proj_disease_ct, Var1, into = c("proj", "disease", "ct"), sep = "_")

col_factors <- c("proj", 
                 "disease",
                 "ct"
)
num_proj_disease_ct[col_factors] <- lapply(num_proj_disease_ct[col_factors], as.factor)  
names(num_proj_disease_ct)[names(num_proj_disease_ct) == 'Freq'] <- "count"

normal_disco <- subset(disco_filt, subset = disease == "Normal")
ad_disco <- subset(disco_filt, subset = disease == "Alzheimer's disease")

DimPlot(normal_disco, reduction="umap", group.by = "disease")
DimPlot(ad_disco, reduction="umap", group.by = "disease")

Idents(normal_disco) <- "ct"

normal_all_degs <- FindAllMarkers(normal_disco, 
            logfc.threshold = 0.25,
            min.pct = 0.1)

ad_all_degs <- FindAllMarkers(normal_disco, 
                                  logfc.threshold = 0.25,
                                  min.pct = 0.1)

write.csv(normal_all_degs, file = paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/", 
                                  "/20220817_DEGs/outputs/normal_all_DEGs.csv"),
                                  row.names = TRUE)

ad_remove <- list.dirs("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/Alzheimer's disease/01A_only_1_project", recursive=FALSE, full.names = FALSE)

# session info
sessionInfo()