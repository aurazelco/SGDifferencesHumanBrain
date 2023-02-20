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
source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Comparison/DEGs/Compare_DEGs_func.R")

# sets the directories where to find the DEG csv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/outputs/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/"

# set the main directory where to save the generated plots - sub-directories are created (if they do not already exist) within the plotting functions
main_comparison <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/"

# Vectors to save the different sub-groups of DISCO and UCSC
sub_disease <- list.dirs(main_DISCO, full.names = F, recursive = F)
# the first folder "exta_files" is excluded
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Threshold to filter the DEGs
pvalue_thresh <- 0.05
FC_thresh <- 1.2

# Import all the CSVs from the different ages/conditions - slightly different file tree structure requires a different approach for UCSC
disco_projs <- c("GSE157827", "GSE174367", "PRJNA544731")
disco <- ImportDataset(main_DISCO, sub_disease, individual_projs = disco_projs, pval = pvalue_thresh, FC = FC_thresh)
UCSC <- ImportDataset(main_UCSC, sub_UCSC, UCSC_flag = "yes")

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
              "Mesenchymal" = "Mesenchymal",      
              "Neuroepithelial" =     "Neuroepithelial",
              "Neuronal" = "Neurons",            
              "Other"    = "Other",                
              "Radial Glial"     = "Radial Glia",       
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
condition_order <- c("Eze_Nowakowski_integrated_2nd_trimester",
               "Velmeshev_2022_2nd_trimester",           
               "Velmeshev_2022_3rd_trimester", 
               "Velmeshev_2022_0_1_years",                
               "Velmeshev_2022_1_2_years",            
               "Velmeshev_2022_2_4_years",  
               "Velmeshev_2022_10_20_years",      
               "Velmeshev_2022_Adult",
               "Normal_GSE157827",              
               "Normal_GSE174367",               
               "Normal_PRJNA544731", 
               "Alzheimer's disease_GSE157827",
               "Alzheimer's disease_GSE174367",
               "Multiple Sclerosis_PRJNA544731" 
)

# Generates a df with all DEGs
sexes <- CreateSexDf(c(UCSC[[1]][-1], disco[[1]]), unified_annotation)

# Heatmaps of presence of genes (yes/no) across all ages, for each ct
PlotCts(main_comparison, sexes, condition_order)

# Count of how genes are shared among ages, for each ct and sex
gene_counts <- CreateCountDfs(sexes)
PlotNumSharedGenes(main_comparison, gene_counts)
SaveSharedGenes(main_comparison, gene_counts, 0.75, 1)

# Count number of DEGs per ct across ages, for each ct and sex, and also create one faceted figure
num_deg <- NumDEGsAcrossConditions(sexes, condition_order[-1])
PlotNumDEGs(main_comparison, num_deg)

brewer_palette <- c(colorRampPalette(c("white", "#228B22"))(8), "#FF0000")
custom_palette <- c(
                    "Velmeshev_2022_2nd_trimester"=brewer_palette[1],           
                    "Velmeshev_2022_3rd_trimester"=brewer_palette[2], 
                    "Velmeshev_2022_0_1_years"=brewer_palette[3],                
                    "Velmeshev_2022_1_2_years"=brewer_palette[4],            
                    "Velmeshev_2022_2_4_years"=brewer_palette[5],  
                    "Velmeshev_2022_10_20_years"=brewer_palette[6],      
                    "Velmeshev_2022_Adult"=brewer_palette[7],
                    "Normal_GSE157827"=brewer_palette[8],              
                    "Normal_GSE174367"=brewer_palette[8],               
                    "Normal_PRJNA544731"=brewer_palette[8], 
                    "Alzheimer's disease_GSE157827"=brewer_palette[9],
                    "Alzheimer's disease_GSE174367"=brewer_palette[9],
                    "Multiple Sclerosis_PRJNA544731"=brewer_palette[9] 
)

PlotNumDEGsFacets(main_comparison, num_deg, custom_palette)

# Plot the number of overlapping genes between one condition and all others, divided by ct and sex
PlotDEGsOverlap(main_comparison, sexes, condition_order)
# Plots the same but as a heatmap
PlotDEGsOverlapHmp(main_comparison, sexes, condition_order)

# Plots to compare 2nd trim
#trim_2nd <- CreateConditionDf(c(UCSC[[1]], disco[[1]]), unified_annotation, condition_order[1:2])
#PlotAcrossConditions(main_comparison, trim_2nd, "trimester_2nd")

# Comparison of normal DISCO DEGs
normal_disco <- NormDf(disco[[1]][c("Normal_GSE157827", "Normal_GSE174367", "Normal_PRJNA544731")], unified_annotation)

og_cts <- CalcCommonGenes(normal_disco, "ct")
common_cts <- CalcCommonGenes(normal_disco, "common_annot")

PlotCommonGenes(main_comparison, og_cts, "Normal_DISCO", "Original_annotation")
PlotCommonGenes(main_comparison, common_cts, "Normal_DISCO", "Unified_annotation")

# Glucocorticoids Binding sites
gbs <- readxl::read_xlsx("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/Polman_2012_GBS.xlsx", skip=9)
gbs_genes <- toupper(gbs$`nearest gene`)

for (i in names(sexes)) {
  for (sex in c("F", "M")) {
    print(paste(i, sex, sep = " "))
    print(length(intersect(gbs_genes, sexes[[i]][which(sexes[[i]]$sex==sex), "gene_id"])))
    print(nrow(sexes[[i]]))
    if (toupper("Nr3c1") %in% sexes[[i]][which(sexes[[i]]$sex==sex), "gene_id"]) {
      print(paste0(toupper("Nr3c1"), " found in ", i, " ", sex))
    } else if (toupper("Nr3c2") %in% sexes[[i]][which(sexes[[i]]$sex==sex), "gene_id"]) {
      print(paste0(toupper("Nr3c2"), " found in ", i, " ", sex))
    } 
  }
}

all_genes <- do.call(rbind, sexes)
all_genes$ct <- gsub("\\..*", "", rownames(all_genes))
all_genes$presence <- str_replace_all(all_genes$presence, c("yes"="Yes", "no"="No"))

grs <- complete(all_genes[which(all_genes$gene_id %in% c( "NR3C1", "NR3C2")),], gene_id, condition,sex,ct )


plot_path <- paste0(main_comparison, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "GR_MR.pdf"), width = 15)
print(
  ggplot(grs, aes(factor(condition, condition_order[which(condition_order %in% unique(condition))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~ ct , scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
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

# Mitochondrial genes

mit_genes_ids <- unique(all_genes$gene_id[which(grepl("^MT-", all_genes$gene_id))])
mit_genes_ids <- c("XIST", mit_genes_ids)

mit_genes <- all_genes[which(all_genes$gene_id %in% mit_genes_ids), ]
mit_gene_count <- as.data.frame(table(mit_genes[which(mit_genes$presence=="Yes"), "gene_id"]))
mit_gene_count <- mit_gene_count[order(mit_gene_count$Freq, decreasing = T),]
mit_genes$gene_id <- factor(mit_genes$gene_id, unique(mit_gene_count$Var1))
mit_genes <- complete(mit_genes, gene_id, condition,sex,ct)

plot_path <- paste0(main_comparison, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "MT_genes.pdf"), width = 15)
print(
  ggplot(mit_genes, 
         aes(factor(condition, condition_order[which(condition_order %in% unique(condition))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~ ct , scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
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


# X-escaping genes

x_escapees <- read.table("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Comparison/escape_Xchr.txt", sep="\t", skip = 2)

plot_path <- paste0(main_comparison, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "X_escaping_genes.pdf"), width = 15)
print(
  ggplot(complete(all_genes[which(all_genes$gene_id %in% x_escapees$V2), ], gene_id, condition,sex,ct), 
         aes(factor(condition, condition_order[which(condition_order %in% unique(condition))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~ ct , scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
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
m_count_F <- as.data.frame(table(all_genes[which(all_genes$presence=="No" & all_genes$sex=="M"), "gene_id"]))
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

f_count_M <- as.data.frame(table(all_genes[which(all_genes$presence=="No" & all_genes$sex=="F"), "gene_id"]))
f_count_M$sex <- rep("F", nrow(f_count_M))
m_count_M <- as.data.frame(table(all_genes[which(all_genes$presence=="Yes" & all_genes$sex=="M"), "gene_id"]))
m_count_M$sex <- rep("M", nrow(m_count_M))
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
most_diff_genes$gene_id <- factor(most_diff_genes$gene_id, rev(c(top10F_u, top10M_u)))

plot_path <- paste0(main_comparison, "Hmp_Presence_Ind_DEGs/")
dir.create(plot_path, recursive = T, showWarnings = F)
pdf(paste0(plot_path, "top_20_most_diff_genes.pdf"), width = 15, height = 15)
print(
  ggplot(most_diff_genes,
         aes(factor(condition, condition_order[which(condition_order %in% unique(condition))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~  ct, scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
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