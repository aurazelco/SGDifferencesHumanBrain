# Author: Aura Zelco
# Brief description:
    # This script is used for comparing the DEGs from the DEG analysis across multiple datasets (different ages/disease conditions)
# Brief procedure:
    # 1. Reads all DEG csv files from all the different datasets (in this case 2 - DISCO and UCSC)
    # 2. Manually combines the annotations to be able to compare at a general level the different celltypes
    # 3. Plots presence heatmaps (yes/no, not the expression) across all ages, for each celltype
    # 4. Plots how many genes are found in all age groups, in all but one, etc
    # 5. Plots the total number of DEGs per ct across all conditions 
    # 6. Plots the number of overlapping genes between a specific condition and all others, divided by ct and sex
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("scripts/Integration/DEGs/DEGs_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "data/DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "data/UCSC/DEGs_adjust_pval/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_int_path <- "data/Integration/"

# Vectors to save the different sub-groups of DISCO and UCSC
sub_disco <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
# the first folder "exta_files" is excluded
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco <- ImportDatasets(main_DISCO, sub_disco, UCSC_flag = "no", individual_projs = T)
names(disco[[1]]) <- str_replace_all(names(disco[[1]]), "Normal", "Healthy")
UCSC <- ImportDatasets(main_UCSC, sub_UCSC, UCSC_flag = "yes", individual_projs = F)

# disco[[2]] and UCSC[[2]] can be used to manually create unified_annotation, as done below
disco[[2]]
UCSC[[2]]

# manually decided how to combine the sub-celltypes
unified_annotation <- c("CXCL14 IN" = "Interneurons",
                        "EC" = "Endothelial cells",
                        "fibrous astrocyte"  = "Astrocytes",
                        "L2_3 EN" = "Excitatory neurons", 
                        "L4 EN" = "Excitatory neurons",
                        "L5 EN" = "Excitatory neurons",
                        "L5_6 EN" = "Excitatory neurons",
                        "L5b EN" = "Excitatory neurons",
                        "L6 EN" = "Excitatory neurons",                
                        "microglia" = "Microglia", 
                        "Oligodendrocyte" =  "Oligodendrocytes",      
                        "OPC" = "OPCs",                  
                        "PLCH1 L4_5 EN" = "Excitatory neurons", 
                        "protoplasmic astrocyte" = "Astrocytes",
                        "PVALB IN"  = "Interneurons",            
                        "pyramidal neuron"  = "Excitatory neurons",
                        "SST IN" = "Interneurons",   
                        "SV2C IN"  = "Interneurons",   
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
                        "Astrocytes" = "Astrocytes",        
                        "Excitatory neurons"  = "Excitatory neurons",
                        "Interneurons"   = "Interneurons",     
                        "Microglia"  = "Microglia",         
                        "Oligodendrocytes" = "Oligodendrocytes",
                        "OPCs" = "OPCs",            
                        "Unknown" = "Unknown",           
                        "Vascular cells" = "Vascular cells",     
                        "Dorsal progenitors"  = "Dorsal progenitors" ,   
                        "Ventral progenitors" = "Ventral progenitors")
names(unified_annotation) <- tolower(names(unified_annotation))

# defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
groups_order <- c(
               "Velmeshev_2022_2nd_trimester",           
               "Velmeshev_2022_3rd_trimester", 
               "Velmeshev_2022_0_1_years",                
               "Velmeshev_2022_1_2_years",            
               "Velmeshev_2022_2_4_years",  
               "Velmeshev_2022_10_20_years",      
               "Velmeshev_2022_Adult",
               "Healthy_GSE157827",              
               "Healthy_GSE174367",               
               "Healthy_PRJNA544731", 
               "Alzheimer's disease_GSE157827",
               "Alzheimer's disease_GSE174367",
               "Multiple Sclerosis_PRJNA544731" 
)

# Generates a df with all DEGs
sexes <- CreateSexDf(c(UCSC[[1]], disco[[1]]), unified_annotation)

# Heatmaps of presence of genes (yes/no) across all ages, for each ct
PlotCts(main_int_path, sexes, groups_order)

# Count of how genes are shared among ages, for each ct and sex
gene_counts <- CreateCountDfs(sexes)
PlotNumSharedGenes(main_int_path, gene_counts)
SaveSharedGenes(main_int_path, gene_counts, 0.75, 1)

PlotDEGsOverlap(main_int_path, sexes, groups_order)

# Count number of DEGs per ct across ages, for each ct and sex, and also create one faceted figure
num_deg <- NumDEGsAcrossConditions(sexes, groups_order)
PlotNumDEGs(main_int_path, num_deg)

brewer_palette <- c(colorRampPalette(c("#DFEEDF", "#1A451E"))(8), "#FF0000")
#brewer_palette <- c(colorRampPalette(c("#44fc57", "#114015"))(8), "#FF0000")
custom_palette <- c(
                    "Velmeshev_2022_2nd_trimester"=brewer_palette[1],           
                    "Velmeshev_2022_3rd_trimester"=brewer_palette[2], 
                    "Velmeshev_2022_0_1_years"=brewer_palette[3],                
                    "Velmeshev_2022_1_2_years"=brewer_palette[4],            
                    "Velmeshev_2022_2_4_years"=brewer_palette[5],  
                    "Velmeshev_2022_10_20_years"=brewer_palette[6],      
                    "Velmeshev_2022_Adult"=brewer_palette[7],
                    "Healthy_GSE157827"=brewer_palette[8],              
                    "Healthy_GSE174367"=brewer_palette[8],               
                    "Healthy_PRJNA544731"=brewer_palette[8], 
                    "Alzheimer's disease_GSE157827"=brewer_palette[9],
                    "Alzheimer's disease_GSE174367"=brewer_palette[9],
                    "Multiple Sclerosis_PRJNA544731"=brewer_palette[9] 
)

PlotNumDEGsFacets(main_int_path, num_deg, custom_palette)

# Plot dotplot for figure 1

cts_order <- c(
  "Excitatory neurons",  
  "Interneurons",        
  "Astrocytes",          
  "Microglia" ,          
  "Oligodendrocytes",    
  "OPCs",   
  "Endothelial cells", 
  "Vascular cells",      
  "Dorsal progenitors",  
  "Ventral progenitors",
  "Unknown"   
)

num_deg$ct <- factor(num_deg$ct, levels = cts_order[which(cts_order %in% unique(num_deg$ct))])

library(gridExtra)
library(gtable)
library(grid)

# Create and extract color legend
p1 <- ggplot(num_deg, aes(condition, ct, fill=condition)) +
  geom_point(color="black", shape=21) +
  scale_fill_manual(values = custom_palette) +
  theme_bw() +
  labs(x="Developmental stages", 
       y="Cell types", 
       fill="Developmental stages") +
  theme(legend.position = "bottom",
        legend.text = element_text(size=14, colour = "black"),
        legend.title = element_text(size=16, face="bold", colour = "black"),) +
  guides(fill = guide_legend(ncol = 3))
leg_color <- gtable_filter(ggplot_gtable(ggplot_build(p1)), "guide-box") 

# Create and extract size legend
p2 <- ggplot(num_deg, aes(condition, ct)) +
  geom_point(aes(size=count_degs)) +
  scale_color_manual(values = custom_palette) +
  scale_y_discrete(limits=rev) +
  theme_bw() +
  labs(x="Developmental stages", 
       y="Cell types", 
       size="Number of DEGs") +
  theme(legend.position = "bottom",
        legend.text = element_text(size=14, colour = "black"),
        legend.title = element_text(size=16, face="bold", colour = "black"),)
leg_size <- gtable_filter(ggplot_gtable(ggplot_build(p2)), "guide-box")

# Create actual plot without any legend
new_fig1 <- ggplot(num_deg, 
       aes(condition, ct, size=count_degs, fill=condition)) + 
  geom_point(color="black", shape=21) +
  scale_fill_manual(values = custom_palette) +
  scale_y_discrete(limits=rev) +
  labs(x="Developmental stages", 
       y="Cell types") +
  facet_wrap(~ sex,drop = T) + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = NA, color = "black"), 
    panel.spacing.x = unit(0.5, "lines"),
    strip.text = element_text(size=14, face="bold", colour = "black"),
    axis.line = element_line(colour = "black"),
    axis.title.y = element_text(size=16, face="bold", colour = "black"),
    axis.text.y = element_text(size=14, colour = "black", vjust = 0.7, hjust=0.5),
    axis.title.x = element_text(size=16, face="bold", colour = "black"),
    axis.text.x = element_blank(), 
    legend.position = "none"
  ) 

# Arrange the three elements in a new plot
fig1_v2 <- ggarrange(new_fig1, leg_size, leg_color, nrow = 3,
          heights = c(1, 0.2, 0.2))

pdf("data/Integration/Num_Total_DEGs_across_conditions/Faceted_tot_DEGs_v2.pdf",
    width = 14, height = 12)
print(fig1_v2)
dev.off()


# Plot the number of overlapping genes between one condition and all others, divided by ct and sex
PlotDEGsOverlap(main_int_path, sexes, groups_order)
# Plots the same but as a heatmap
PlotDEGsOverlapHmp(main_int_path, sexes, groups_order)

# Comparison of Healthy DISCO DEGs
Healthy_disco <- NormDf(disco[[1]][c("Healthy_GSE157827", "Healthy_GSE174367", "Healthy_PRJNA544731")], unified_annotation)

og_cts <- CalcCommonGenes(Healthy_disco, "ct")
og_cts$ct <- str_to_title(og_cts$ct)
common_cts <- CalcCommonGenes(Healthy_disco, "common_annot")

PlotCommonGenes(main_int_path, og_cts, "Healthy_DISCO", "Original_annotation")
PlotCommonGenes(main_int_path, common_cts, "Healthy_DISCO", "Unified_annotation")

og_cts$annot_type <- rep("Original annotation", nrow(og_cts))
common_cts$annot_type <- rep("Unified annotation", nrow(common_cts))

all_annot <- rbind(og_cts, common_cts)

PlotCommonGenes(main_int_path, all_annot, "Healthy_DISCO", "both")

# Mitochondrial genes

mit_genes_ids <- unique(all_genes$gene_id[which(grepl("^MT-", all_genes$gene_id))])
timtom_ids <- c(unique(all_genes$gene_id[which(grepl("^TIMM", all_genes$gene_id))]),
                unique(all_genes$gene_id[which(grepl("^TOMM", all_genes$gene_id))]))
tca_genes <- c("ACO2", "CS", "FH", "MDH1", "OGDH", "PDHA1", "PDHA2", "SDHC", "SUCLG1")
# https://maayanlab.cloud/Harmonizome/gene_set/TCA+cycle/PANTHER+Pathways
mit_genes_ids <- c("XIST", timtom_ids, tca_genes, mit_genes_ids)

mit_genes <- all_genes[which(all_genes$gene_id %in% mit_genes_ids), ]
mit_gene_count <- as.data.frame(table(mit_genes[which(mit_genes$presence=="Yes"), "gene_id"]))
mit_gene_count <- mit_gene_count[order(mit_gene_count$Freq, decreasing = T),]
mit_genes$gene_id <- factor(mit_genes$gene_id, unique(mit_gene_count$Var1))
mit_genes <- complete(mit_genes, gene_id, condition,sex,ct)

mit_order <- c(
  "XIST", "TIMM17A", "TIMM23B", "TOMM7", "TOMM20", "MT-CYB", "MT-ND1",  
  "MT-ND2",  "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-CO1", "MT-CO2",  
  "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-RNR1", "MT-RNR2", "ACO2", "CS", "FH", 
  "MDH1", "OGDH", "PDHA1", "PDHA2", "SDHC", "SUCLG1")

plot_path <- paste0(main_int_path, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "MT_genes.pdf"), width = 9, height = 14)
print(
  ggplot(mit_genes, 
         aes(factor(gene_id, rev(unique(mit_order))), factor(condition, rev(groups_order[which(groups_order %in% unique(condition))])),  fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(ct ~ sex , scales = "free", space = "free") +
    labs(y="Datasets", x="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0.5, "lines"),
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  
  
)
dev.off()


dis_mit_genes <- readxl::read_xlsx("Integration/Steinmetz_suppl_2022.xlsx", sheet=1, skip=1)

dis_genes <- unique(dis_mit_genes$`gene name`)


dis_mit <- all_genes[which(all_genes$gene_id %in% dis_genes), ]
dis_mit_gene_count <- as.data.frame(table(dis_mit[which(dis_mit$presence=="Yes"), "gene_id"]))
dis_mit_gene_count <- dis_mit_gene_count[order(dis_mit_gene_count$Freq, decreasing = T),]
dis_mit$gene_id <- factor(dis_mit$gene_id, unique(dis_mit_gene_count$Var1))
dis_mit <- complete(dis_mit, gene_id, condition,sex,ct)

plot_path <- paste0(main_int_path, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "MT_genes_Steinmetz.pdf"), width = 9, height = 14)
print(
ggplot(dis_mit, 
       aes(gene_id,factor(condition, groups_order[which(groups_order %in% unique(condition))]),  fill=presence)) +
  geom_tile() +
  scale_fill_manual(values = c("Yes"="#F8766D",
                               "No"="#00BFC4"),
                    na.value = "#00BFC4",
                    guide = guide_legend(reverse = TRUE)) +
  facet_grid(ct ~ sex , scales = "free", space = "free") +
  labs(y="Groups", x="Genes", fill="Genes found") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.spacing.x=unit(0, "lines"),
        strip.text.x = element_text(size=12, face="bold", colour = "black"),
        strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()





# X-escaping genes

x_escapees <- read.table("Integration/escape_Xchr.txt", sep="\t", skip = 2)

plot_path <- paste0(main_int_path, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "X_escaping_genes.pdf"), width = 9, height = 14)
print(
  ggplot(complete(all_genes[which(all_genes$gene_id %in% x_escapees$V2), ], gene_id, condition,sex,ct), 
         aes(gene_id, factor(condition, rev(groups_order[which(groups_order %in% unique(condition))])),  fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(ct ~ sex, scales = "free", space = "free") +
    labs(y="Datasets", x="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  
  
)
dev.off()

# Genes with high difference between M and F

f_count_F <- as.data.frame(table(all_genes[which(all_genes$presence=="Yes" & all_genes$sex=="F"), "gene_id"]))
f_count_F$sex <- rep("F", nrow(f_count_F))
m_count_F <- as.data.frame(table(all_genes[which((all_genes$gene_id %in% unique(f_count_F$Var1)) & all_genes$sex=="M" & all_genes$presence=="Yes"), "gene_id"]))
m_count_F$sex <- rep("M", nrow(m_count_F))
sex_count_F <- rbind(f_count_F, m_count_F)
colnames(sex_count_F) <- c("gene_id", "count", "sex")
sex_count_F <- complete(sex_count_F, gene_id, sex)
sex_count_F[which(is.na(sex_count_F$count)), "count"] <- 0

abs_diff_F <- data.frame("gene_id"=unique(sex_count_F$gene_id))
abs_diff_F$sex_diff <- rep(NA, nrow(abs_diff_F))
for (i in abs_diff_F$gene_id) {
  abs_diff_F[which(abs_diff_F$gene_id==i), "sex_diff"] <- sex_count_F[which(sex_count_F$sex=="F" & sex_count_F$gene_id==i), "count"] - sex_count_F[which(sex_count_F$sex=="M" & sex_count_F$gene_id==i), "count"]
}
abs_diff_F <- abs_diff_F[which(abs_diff_F$sex_diff > 0 ),]
abs_diff_F <- abs_diff_F[order(abs_diff_F$sex_diff, decreasing = T), ]
abs_diff_F$gene_id <- factor(abs_diff_F$gene_id, unique(abs_diff_F$gene_id))


m_count_M <- as.data.frame(table(all_genes[which(all_genes$presence=="Yes" & all_genes$sex=="M"), "gene_id"]))
m_count_M$sex <- rep("M", nrow(m_count_M))
f_count_M <- as.data.frame(table(all_genes[which((all_genes$gene_id %in% unique(m_count_M$Var1)) & all_genes$sex=="F" & all_genes$presence=="Yes"), "gene_id"]))
f_count_M$sex <- rep("F", nrow(f_count_M))
sex_count_M <- rbind(f_count_M, m_count_M)
colnames(sex_count_M) <- c("gene_id", "count", "sex")
sex_count_M <- complete(sex_count_M, gene_id, sex)
sex_count_M[which(is.na(sex_count_M$count)), "count"] <- 0

abs_diff_M <- data.frame("gene_id"=unique(sex_count_M$gene_id))
abs_diff_M$sex_diff <- rep(NA, nrow(abs_diff_M))
for (i in abs_diff_M$gene_id) {
  abs_diff_M[which(abs_diff_M$gene_id==i), "sex_diff"] <- sex_count_M[which(sex_count_M$sex=="M" & sex_count_M$gene_id==i), "count"] - sex_count_M[which(sex_count_M$sex=="F" & sex_count_M$gene_id==i), "count"]
}
abs_diff_M <- abs_diff_M[which(abs_diff_M$sex_diff > 0 ),]
abs_diff_M <- abs_diff_M[order(abs_diff_M$sex_diff, decreasing = T), ]
abs_diff_M$gene_id <- factor(abs_diff_M$gene_id, unique(abs_diff_M$gene_id))

top10F <- as.character(abs_diff_F[1:10, "gene_id"])
top10M <- as.character(abs_diff_M[1:10, "gene_id"])

common_genes <- intersect(top10F, top10M)

`%!in%` <- Negate(`%in%`)

if (length(common_genes) > 0 ) {
  top10F_u <- top10F[which(top10F %!in% common_genes)]
  top10M_u <- top10M[which(top10M %!in% common_genes)]
  while (length(top10F_u) < 10 & length(top10M_u) < 10) {
    print(length(top10F_u))
    print(length(top10M_u))
    for (i  in common_genes) {
      if (abs_diff_F[which(abs_diff_F$gene_id==i), "sex_diff"] > abs_diff_M[which(abs_diff_M$gene_id==i), "sex_diff"]) {
        top10F_u <- c(top10F_u, i)
      } else {
        top10M_u <- c(top10M_u, i)
      }
    }
  }
} else {
  top10F_u <- top10F
  top10M_u <- top10M
}

if (length(top10F_u) < 10) {
  missing_genes <- as.character(abs_diff_F$gene_id[11:(11+10 - length(top10F_u) - 1)])
  top10F_u <- c(top10F_u, missing_genes)
}
if (length(top10M_u) < 10) {
  missing_genes <- as.character(abs_diff_M$gene_id[11:(11+10 - length(top10M_u) - 1)])
  top10M_u <- c(top10M_u, missing_genes)
}


most_diff_genes <- c(top10F_u, top10M_u)
most_diff_genes <- complete(all_genes[which(all_genes$gene_id %in% most_diff_genes), ], gene_id, condition,sex,ct)
most_diff_genes$gene_id <- factor(most_diff_genes$gene_id, c(top10F_u, top10M_u))

plot_path <- paste0(main_int_path, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "top_20_most_diff_genes.pdf"), width = 9, height = 15)
print(
  ggplot(most_diff_genes,
         aes(gene_id, factor(condition, rev(groups_order[which(groups_order %in% unique(condition))])), fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(ct ~  sex, scales = "free", space = "free") +
    labs(y="Datasets", x="Genes", fill="SG-biased genes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0.5, "lines"),
          strip.text.x = element_text(size=12, face="bold", colour = "black"),
          strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  
)
dev.off()



library(tidyverse)
library(purrr)
library(RColorBrewer)


disco_folders <- sapply(sub_disco, function(x) {
  deg_folders <- list.dirs(paste0(main_DISCO, x), full.names = T, recursive = T)
  deg_folders <- deg_folders[which(grepl("01B_num_DEGs", deg_folders))]
})
  
ucsc_folders <- sapply(sub_UCSC, function(x) {
  deg_folders <- list.dirs(paste0(main_UCSC, x), full.names = T, recursive = T)
  deg_folders <- deg_folders[which(grepl("01B_num_DEGs", deg_folders))]
})

deg_folders <- unlist(c(ucsc_folders, disco_folders), use.names = F)

deg_list <- lapply(deg_folders, function(x) {
  ImportDE(x)
})
names(deg_list) <- deg_folders

deg_list <- deg_list[lapply(deg_list,length)>0]


names(deg_list) <- str_remove_all(names(deg_list), paste(c("data/DISCOv1.0/DEGs_proj_adjust_pval/",
                                                         "data/UCSC/DEGs_adjust_pval/",
                                                         "/outputs/01B_num_DEGs"), 
                                                         collapse = "|"))
deg_names <- str_replace_all(names(deg_list), 
                             c("GSE174367/" = "GSE174367_", 
                               "GSE157827/"="GSE157827_", 
                               "PRJNA544731/"="PRJNA544731_", 
                               "Normal"="Healthy"))

deg_list <- lapply(1:length(deg_list), function(x) {
  tmp <- lapply(names(deg_list[[x]]), function(y) {
    deg_list[[x]][[y]] %>%
      rownames_to_column(var = "gene") %>%
      select(-index)
  })
  names(tmp) <- str_remove_all(names(deg_list[[x]]), "_filt")
  new_list <- bind_rows(tmp, .id = "sex")
  return(new_list)
})
names(deg_list) <- deg_names

annot_df <- tibble(names(unified_annotation), unified_annotation) %>%
            rename("ct" = `names(unified_annotation)`)

full_df <- bind_rows(deg_list, .id = "groups_id") %>%
            separate(groups_id, into = c("group", "ct"), sep = "/") %>%
            mutate(ct = tolower(ct)) %>%
            left_join(., annot_df, by = "ct") %>%
            select(-ct) %>%
            rename("ct" = unified_annotation) %>%
            select(group, ct, everything())

mit_genes_ids <- full_df %>% 
                  filter(grepl("^MT-|^TIMM|^TOMM", gene)) %>% 
                  pull(gene) %>% 
                  unique()
tca_genes <- c("ACO2", "CS", "FH", "MDH1", "OGDH", "PDHA1", "PDHA2", "SDHC", "SUCLG1")
# https://maayanlab.cloud/Harmonizome/gene_set/TCA+cycle/PANTHER+Pathways
mit_genes_ids <- c("XIST", mit_genes_ids, tca_genes)
mit_order <- c(
  "XIST", "TIMM17A", "TIMM23B", "TOMM7", "TOMM20", "MT-CYB", "MT-ND1",  
  "MT-ND2",  "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-CO1", "MT-CO2",  
  "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-RNR1", "MT-RNR2", "ACO2", "CS", "FH", 
  "MDH1", "OGDH", "PDHA1", "PDHA2", "SDHC", "SUCLG1")



check_duplicates <- full_df %>% 
                      filter(gene %in% mit_genes_ids) %>%
                      group_by(group, sex, gene) %>% 
                      count() %>%
                      filter(n > 1)



mt_df <- full_df %>%
          filter(gene %in% mit_genes_ids) %>%
          mutate(gene = factor(gene, levels = mit_order),
                 group = factor(
                           str_replace_all(group, c(
                             "GSE157827_Alzheimer's disease" = "Alzheimer's disease_GSE157827", 
                             "GSE157827_Healthy" = "Healthy_GSE157827",
                             "GSE174367_Alzheimer's disease" = "Alzheimer's disease_GSE174367",
                             "GSE174367_Healthy" = "Healthy_GSE174367", 
                             "PRJNA544731_Multiple Sclerosis" = "Multiple Sclerosis_PRJNA544731",
                             "PRJNA544731_Healthy" = "Healthy_PRJNA544731" 
                           )), 
                          levels = groups_order),
                 pval_sign = case_when(
                               p_val_adj > 0.05 ~ "NS",
                               p_val_adj > 0.01 & p_val_adj <= 0.05 ~ "*",
                               p_val_adj > 0.001 & p_val_adj <= 0.01 ~ "**",
                               p_val_adj > 0.0001 & p_val_adj <= 0.001 ~ "***",
                              p_val_adj <= 0.0001 ~ "****"
                                       ))




brewer_palette <- brewer.pal(4,"Purples")

p <- ggplot(mt_df, aes(gene, group, fill=pval_sign, size=avg_log2FC)) +
  geom_point(color="black", shape=21) +
  facet_grid(ct ~ sex, scales = "free_y", space = "free") +
  scale_fill_manual(values = c("NS"="white", 
                               "*"=brewer_palette[1],
                               "**"=brewer_palette[2],
                               "***"=brewer_palette[3],
                               "****"=brewer_palette[4]),
                    na.value = "gray") +
  labs(x="Genes of interest", y="Groups", size="Average Log2FC", fill="Adjusted p-values") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        panel.spacing.x=unit(0.5, "lines"),
        strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
        strip.text.x = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
        axis.text.y = element_text(size=8, colour = "black"),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.box = "vertical",
        legend.title = element_text(size=12, face="bold", colour = "black"))



pdf(paste0(main_int_path, "Hmp_Presence_Ind_DEGs/MT_genes_v2.pdf"), width = 10, height = 14)
print(p)
dev.off()


library(biomaRt)

mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "useast")
Annot_idf <- getBM(attributes = c("hgnc_symbol",
                                  "chromosome_name",
                                  "description"),
                   filters = c("hgnc_symbol") ,
                   values=list(full_df %>% pull(gene) %>% unique()),
                   mart = mart)
mt_annot <- Annot_idf %>%
              filter(chromosome_name == "MT")



# Ribosomial genes

library(dplyr)

rb_genes_ids <- unique(all_genes$gene_id[which(grepl("^RPL|^RPS|^MRP", all_genes$gene_id))])
rb_genes <- all_genes[which(all_genes$gene_id %in% rb_genes_ids), ]
rb_gene_count <- as.data.frame(table(rb_genes[which(rb_genes$presence=="Yes"), "gene_id"]))
rb_gene_count <- rb_gene_count[order(rb_gene_count$Freq, decreasing = T),]
rb_genes$gene_id <- factor(rb_genes$gene_id, unique(rb_gene_count$Var1))
rb_genes <- complete(rb_genes, gene_id, condition,sex,ct)

rb_genes <- rb_genes %>%
              mutate("ribo_group" = case_when(grepl("^RPL", gene_id) ~ "RPL",
                                              grepl("^RPS", gene_id) ~ "RPS",
                                              grepl("^MRP", gene_id) ~ "MRP"))

plot_path <- paste0(main_int_path, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)

lapply(unique(rb_genes$ribo_group), function(x) {
  pdf(paste0(plot_path, "ribo_genes_", x,".pdf"), width = 12, height = 15)
  print(
    ggplot(rb_genes %>% filter(ribo_group == x), 
           aes(gene_id, factor(condition, rev(groups_order[which(groups_order %in% unique(condition))])),  fill=presence)) +
      geom_tile() +
      scale_fill_manual(values = c("Yes"="#F8766D",
                                   "No"="#00BFC4"),
                        na.value = "#00BFC4",
                        guide = guide_legend(reverse = TRUE)) +
      facet_grid(ct ~ sex, scales = "free", space = "free") +
      labs(y="Datasets", x="Genes", fill="Genes found") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            panel.spacing.x=unit(0.5, "lines"),
            strip.text.x = element_text(size=12, face="bold", colour = "black"),
            strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
            plot.title = element_text(size=12, face="bold", colour = "black"),
            axis.line = element_line(colour = "black"),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
    
    
  )
  dev.off()
  
  
})

