library(Seurat)
library(SeuratObject)
library(tidyr)
library(stringr)
library(ggplot2)

rds_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Conde_2022/"


##################################### 1. Filter cts for DEGs analysis steps

conde <- readRDS(paste0(rds_path, "conde_2022.rds"))

# retrieves the number of donors per tissue and sex

conde@meta.data$organ_sex <- paste(conde@meta.data$Organ, conde@meta.data$Sex, sep = "--")

num_samples <- as.data.frame(unique(conde@meta.data$organ_se))
colnames(num_samples) <- c("id")
num_samples <- separate(num_samples, id, into = c("organ", "sex"), sep = "--", remove = F)

num_samples$num_donors <- rep(NA, nrow(num_samples))
for (id in num_samples$id) {
  num_samples[which(num_samples$id==id), "num_donors"] <- length(unique(conde@meta.data[which(conde@meta.data$organ_sex==id), "Donor"]))
}

num_samples <- subset(num_samples, num_donors>2)

orgs_to_keep <- vector()
for (org in unique(num_samples$organ)) {
  if (length(unique(num_samples[which(num_samples$organ==org), "sex"]))==2) {
    orgs_to_keep <- c(orgs_to_keep, org)
  }
}

num_samples <- num_samples[which(num_samples$organ %in% orgs_to_keep), c(2:4)]

write.csv(num_samples, paste0(rds_path, "outputs/num_samples_filt.csv"))

filt_organs <- unique(num_samples$organ)


# extract the number of cells per combination of project, sex, disease condition and celltype
num_sex_organ_ct <- as.data.frame(table(conde@meta.data$sex_organ_ct))
num_sex_organ_ct <- separate(num_sex_organ_ct, Var1, into = c("sex" , "organ", "ct"), sep = "--")
num_sex_organ_ct <- subset(num_sex_organ_ct, organ %in% filt_organs)
col_factors <- c("sex" , "organ", "ct")
num_sex_organ_ct[col_factors] <- lapply(num_sex_organ_ct[col_factors], as.factor)  
names(num_sex_organ_ct)[names(num_sex_organ_ct) == 'Freq'] <- "count"



# filters the above df to check which celltypes contain less than a custom threshold number of cells
# then removes also the counterpart sex if only one of the 2 is less than the threshold
FiltDF <- function(df, min_num_cells) {
  `%!in%` <- Negate(`%in%`)
  df <- droplevels(df)
  incomplete_proj <- vector()
  for (type in levels(df$ct)) {
    for (id in levels(df$organ)) {
      if ((nrow(subset(df, subset = (ct==type & organ==id)))%%2!=0) | (any(subset(df, subset = (ct==type & organ==id))[,"count"] < min_num_cells))) {
        incomplete_proj <- c(incomplete_proj, (paste(id, type, sep="_")))
      }
    }
  }
  df$og <- paste(df$organ, df$ct, sep="_")
  df_filt <- df[df$og %!in% incomplete_proj,]
  return(df_filt)
}

# splits the dfs in three subsets, one per disease condition

min_num_cells <- c(10,50,100)

degs_path <- paste0(rds_path,"DEGs/")
dir.create(degs_path, showWarnings = F, recursive = T)

# saves the cSV outputs for the three thresholds
for (min_cells in min_num_cells) {
  num_filt <- FiltDF(num_sex_organ_ct, min_cells)
  num_filt$idents <- paste(num_filt$sex, num_filt$organ, num_filt$ct, sep="_")
  num_filt$name_subfolders <- str_replace_all(num_filt$ct, "/", "_")
  write.csv(num_filt, file = paste0(degs_path, "final_filt_", min_cells, ".csv"),
            row.names = F)
}

# and the corresponding plots
for (min_cells in min_num_cells) {
  pdf(paste0(degs_path, "filt_counts_", min_cells, ".pdf"), width = 15,  height = 15)
  print(ggplot(num_sex_organ_ct, aes(ct, count, fill=sex)) +
          geom_bar(stat="identity", position = "dodge") + 
          labs(x="", y="Nuclei count", fill="Sex") +
          facet_wrap(~ organ, scales = "free") +
          geom_hline(yintercept = min_cells, linetype="dashed") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.title.x = element_text(size=12, face="bold", colour = "black"),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.5, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom")
  )
  dev.off()
}

# updates the RDS
saveRDS(conde, paste0(rds_path, "conde_2022.rds"))


############ For 02C_Conservation - DEGs

expr_mat_all_cts <- GetAssayData(conde[["RNA"]], slot="data")

rm(conde)
expr_mat_all_cts <- as.data.frame(as.matrix(expr_mat_all_cts))

conde <- readRDS(paste0(rds_path, "conde_2022.rds"))

Idents(conde) <- "sex_organ_ct"



# creates a df containing metadata associated with the cell barcodes -> used later to create the randomly sampled matrices
cell_info_ls <- list()

for (org in orgs_to_keep) {
  org_obj <- subset(conde, Organ==org)
  cell_info <- data.frame()
  for (i in unique(org_obj@meta.data$sex_organ_ct)) {
    print(i)
    cell_id <- WhichCells(org_obj, idents = i)
    og_group <- rep(i, length(cell_id))
    cell_info <- rbind(cell_info, data.frame(cell_id, og_group))
  }
  org_path <- paste0(degs_path, org, "/")
  dir.create(org_path, showWarnings = F, recursive = T)
  cell_info_ls <- append(cell_info_ls, list(cell_info))
  write.csv(cell_info, paste0(org_path, "cell_info.csv"))
}
names(cell_info_ls) <- orgs_to_keep

cell_info_ls <- do.call(rbind, cell_info_ls)

df_list <- list()
df_list_n <- vector()
for (df_id in unique(cell_info_ls$og_group)) {
  og_cells <- c("Genes", cell_info_ls[which(cell_info_ls$og_group==df_id), "cell_id"])
  df_og_df <- expr_mat_all_cts[ , (names(expr_mat_all_cts) %in% og_cells)]
  df_list <- append(df_list, list(df_og_df))
  df_list_n <-  c(df_list_n, df_id)
}
names(df_list) <- df_list_n

rm(cell_info, df_og_df, df_list_n, org_obj)

sub_orgs <- list.dirs(degs_path, recursive=FALSE, full.names = FALSE)
sub_orgs <- sub_orgs[c(1, 3:6)]


dfs_org <- list()
for (org in sub_orgs) {
  dfs_org <- append(dfs_org, list(df_list[names(df_list)[which(grepl(org, names(df_list)))]]))
}
names(dfs_org) <- sub_orgs

FiltOrganDf <- function(df_list) {
  filt_names <- vector()
  df_dis <- list()
  for (k in names(df_list)) {
    if (!is.null(ncol(df_list[[k]]))) {
      df_dis <- append(df_dis, list(rownames(df_list[[k]][which(rowSums(as.matrix(df_list[[k]]))!=0),])))
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

tot_df <- lapply(1:length(names(dfs_org)), function(x) FiltOrganDf(dfs_org[[x]]))
tot_df <- do.call(rbind, tot_df)

tot_df <- separate(tot_df, cts, into=c("sex", "organ", "ct"), sep ="--", remove = FALSE)
names(tot_df)[names(tot_df) == "cts"] <- "og"


write.csv(tot_df, paste0(degs_path, "tot_genes_ct.csv"))

######## DEGs

num_filt <- read.csv(paste0(degs_path, "final_filt_100.csv"))
col_factors <- c("sex", "organ", "ct", "og", "idents")
num_filt[col_factors] <- lapply(num_filt[col_factors], as.factor)  

num_filt$ct <- factor(num_filt$ct, levels=c(
  "Classical monocytes",  
  "Nonclassical monocytes",   
  "Alveolar macrophages",     
  "Intermediate macrophages MNP/T", 
  "Cycling T&NK",             
  "NK_CD16+",     
  "NK_CD56bright_CD16-",      
  "Teffector/EM_CD4",  
  "Trm/em_CD8",
  "Tregs",                    
  "ABCs", 
  "Memory B cells",          
  "Plasma cells",  
  "Plasmablasts",             
  "doublets"           
))

for (org in levels(num_filt$organ)) {
  dir.create(paste0(degs_path, org), showWarnings = F)
}

# Plot number of cells per sex and ct
for (org in levels(num_filt$organ)) {
  pdf(paste0(degs_path, org, "/num_ct_sex.pdf"))
  print(
    ggplot(num_filt[which(num_filt$organ==org), ], aes(ct, count, fill=ct)) +
      geom_bar(stat="identity") +
      facet_wrap(~sex, scales="free") +
      geom_hline(yintercept = 100, linetype="dashed") +
      labs(y="Number of cells", fill="") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom")
  )  
  dev.off()
}

