setwd("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/")
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(data.table)
library(dplyr)
library(matrixStats)
library(stringr)
`%!in%` <- Negate(`%in%`)

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

write.csv(normal_all_degs, file = paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/", 
                                  "/20220817_DEGs/outputs/normal_all_DEGs.csv"),
                                  row.names = TRUE)

Idents(ad_disco) <- "ct"

ad_all_degs <- FindAllMarkers(ad_disco, 
                              logfc.threshold = 0.25,
                              min.pct = 0.1)

write.csv(ad_all_degs, file = paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/", 
                                         "/20220817_DEGs/outputs/ad_all_DEGs.csv"),
          row.names = TRUE)

#ad_remove <- list.dirs("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/Alzheimer's disease/01A_only_1_project", recursive=FALSE, full.names = FALSE)


####### INFO FOR SCENIC PIPELINE

disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

num_proj_sex <- as.data.frame(table(disco_filt$proj_sex_disease))
num_proj_sex <- separate(num_proj_sex, Var1, into = c("proj", "sex", "disease"), sep = "_")
names(num_proj_sex)[names(num_proj_sex) == 'Freq'] <- "count"
write.csv(num_proj_sex, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/num_proj_sex_disease.csv")

for (dis_type in unique(num_proj_sex$disease)) {
  pdf(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/num_proj_sex_", dis_type,".pdf"))
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
write.csv(num_proj_sex_ct, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/num_proj_sex_ct.csv")

for (dis_type in unique(num_proj_sex_ct$disease)) {
  pdf(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/num_proj_sex_", dis_type,"_ct.pdf"),
      width = 15, height = 18.75)
  print(
      ggplot(num_proj_sex_ct[which(num_proj_sex_ct$disease==dis_type), ], aes(proj, count, fill=sex)) +
        geom_bar(stat="identity", position="dodge") +
        labs(x="Project ID", y='Cell Counts', fill="Sex", main = dis_type) +
        facet_wrap(~ct, scales = "free") +
        geom_hline(yintercept = 10, linetype=2) +
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

######## HIGHEST VARIABLE GENES - TEST 2022.10.17

disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

#length(test1$RNA)
#length(test1$RNA[which(test1$RNA>0),])

Idents(disco_filt) <- "proj_sex_disease_ct"


test2 <- AverageExpression(disco_filt, assays="RNA") # calculates average expression for each gene in each ident
test3 <- as.data.frame(test2$RNA) # creates a df

library(matrixStats)
test3$SD <- rowSds(as.matrix(test3)) # calculates SD for each row -> variability in the different groups
test3 <- test3[order(-test3$SD),] # order df in descending order
test3$Genes <- rownames(test3) # extract gene names as factor
test3$Genes <- as.factor(test3$Genes)

# displays SD as descending values
ggplot(test3[1:50,], aes(reorder(Genes, -SD), SD)) +
  geom_point()


# this theoretically retrieves the avg expression of genes within a cluster -> variability among samples
micro1 <- subset(disco_filt, subset = proj_sex_disease_ct == "GSE174367_F_Normal_microglia")

Idents(micro1) <- "sample_id"
test1 <- AverageExpression(micro1, assays="RNA")
test1 <- as.data.frame(test1$RNA) # creates a df

test1$SD <- rowSds(as.matrix(test1)) # calculates SD for each row -> variability in the different groups
test1 <- test1[order(-test1$SD),] # order df in descending order
test1$Genes <- rownames(test1) # extract gene names as factor
test1$Genes <- as.factor(test1$Genes)

# displays SD as descending values
ggplot(test1[1:50,], aes(reorder(Genes, -SD), SD)) +
  geom_point()

test2 <- as.data.frame(GetAssayData(micro1[["RNA"]], slot="data"))
test2$SD <- rowSds(as.matrix(test2))
test2 <- test2[order(-test2$SD),] # order df in descending order
test2$Genes <- rownames(test2) # extract gene names as factor
test2$Genes <- as.factor(test2$Genes)

# displays SD as descending values
ggplot(test2[1:50,], aes(reorder(Genes, -SD), SD)) +
  geom_point()


# main diff is that AverageExpression is in log-space, GetAssayData is not -> use GetAssayData

sub_test2 <- test2[1:2000,]
drops <- c("Genes", "SD")
sub_test2 <- sub_test2[ , !(names(sub_test2) %in% drops)]
sub_test2 <- cbind("Genes" = rownames(sub_test2), sub_test2)
rownames(sub_test2) <- NULL


sub_test2 <- rbind(colnames(sub_test2), sub_test2)

write.table(sub_test2, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/test_micro.tsv", 
            sep="\t", col.names = FALSE)
write.table(sub_test2, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/microglia_test/test_micro.tsv", 
            sep="\t", col.names = FALSE)

sub_test3 <- as.data.frame(t(sub_test2))
rownames(sub_test3) <- NULL
sub_test3[1,1] <- NA
write.table(sub_test3, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/microglia_test/test_micro.csv", 
            sep=",", col.names = FALSE, row.names = FALSE)

library(SeuratDisk)
lfile <- SaveLoom(micro1, 
        filename = "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/microglia_test/micro_test.loom")

######## HIGHEST VARIABLE GENES - 2022.10.18

disco_filt <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

#disco_filt@meta.data$proj_sex_disease_ct_sample <- paste(disco_filt@meta.data$proj_sex_disease_ct, disco_filt@meta.data$sample_id, sep="_")
#saveRDS(disco_filt, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/brainV1.0_all_FM_filt.rds")

Idents(disco_filt) <- "proj_sex_disease_ct_sample"

expr_mat_all <- GetAssayData(disco_filt[["RNA"]], slot="data")


Idents(disco_filt) <- "proj_sex_disease_ct"

cell_info <- data.frame()

for (i in unique(disco_filt@meta.data$proj_sex_disease_ct)) {
  #print(i)
  cell_id <- WhichCells(disco_filt, idents = i)
  og_group <- rep(i, length(cell_id))
  cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
}

cell_info <- separate(cell_info, cell_id, into = c("barcode", "sample"), sep = "--", remove = FALSE)

write.csv(cell_info,
          "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/cell_info.csv")

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

saveRDS(expr_mat_all, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/top_2000_SD_expr_matrix.rds")

##### Extract randomly 2k cells and analyze them in SCENIC to see if it can handle a big dataset
random_expr <- sample(expr_mat_all[-2], 2000)
random_expr <- cbind("Genes" = expr_mat_all$Genes, random_expr)
random_expr <- rbind(colnames(random_expr), random_expr)

random_expr_tiny <- random_expr[, 1:1000]

write.table(random_expr, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/random_expr_test/random_expr.csv", 
            sep=",", col.names = FALSE, row.names = FALSE)
write.table(random_expr_tiny, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/SCENIC/random_expr_test/random_expr_tiny.csv", 
            sep=",", col.names = FALSE, row.names = FALSE)


##### Map the samples back to the groups they belong to
#cell_info$cell_id <- paste(cell_info$barcode, cell_info$sample, sep = "--")
#cell_info <- cell_info %>% 
#  relocate(cell_id, .before = barcode)

expr_mat_all <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/top_2000_SD_expr_matrix.rds")



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
  pdf(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/num_filt_", thresh, "_cells.pdf"),
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

expr_mat_all <- readRDS("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/top_2000_SD_expr_matrix.rds")
cell_info <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/cell_info.csv")

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

maindir <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20221018_SCENIC/"

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
  dir.create(paste0(main, dis_type, "/sampled_", num_cells, "_cells"), showWarnings = FALSE)
  lapply(1:length(names(sampled_dfs)), function(i) write.csv(sampled_dfs[[i]], 
                                                              file = paste0(main, dis_type, "/sampled_", num_cells, "_cells/", names(sampled_dfs)[i], ".csv"),
                                                              row.names = FALSE))
  return(sampled_dfs)
}


norm_sampled <- rand_sample(norm, 3, 100, maindir, sub_disease[3])
ad_sampled <- rand_sample(ad, 3, 100, maindir, sub_disease[1])
ms_sampled <- rand_sample(ms, 3, 100, maindir, sub_disease[2])


# session info
sessionInfo()