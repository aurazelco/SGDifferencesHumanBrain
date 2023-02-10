# Author: Aura Zelco
  # Brief description:
    # This script is used for calculating the number of samples used and the XIST expression
  # Brief procedure:
    # 1. reads the RDS and saves the RNA expression for XIST for DISCO and Velmeshev dataset (10-20 years on server, 2nd trimester no XIST)
    # 2. Merges the XIST results into a df and adds the metadata imported at the same time as the RDSs
    # 3. Calculates the number of samples for paper
    # 4. Plots the XIST expression
  # Documentation abbreviations:
    # F and M: females and males
    # ct: celltype
    # df: dataframe
    # ds: dataset

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

library(Seurat)
library(stringr)
library(tidyr)
library(ggplot2)

#---------------------------------------------------------------------------------------------------
####### ON FURU
out_path <- "/Home/ii/auraz/data/UCSC/XIST/"
dir.create(out_path, recursive = T, showWarnings = F)

# Velmeshev 2nd trimester -> No XIST

# Velmeshev 10-20 years
velm_10_20_years <- readRDS("/Home/ii/auraz/data/UCSC/Seurat_UCSC/Velmeshev/Velmeshev_2022_10_20_years.rds")
velm_10_20_years_xist <- GetAssayData(velm_10_20_years[["RNA"]], slot="data")["XIST",]
write.csv(as.data.frame(velm_10_20_years_xist), paste0(out_path, "Velmeshev_2022_10_20_years_XIST.csv"))

#---------------------------------------------------------------------------------------------------
####### LOCAL

# 1. Read local RDS and save XIST info in list, and cell info in another list

ds_paths <- c("Velmeshev_3rd_trimester"= "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/",
             "Velmeshev_0_1_years"="/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/",
             "Velmeshev_1_2_years"= "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/",
             "Velmeshev_2_4_years"= "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/",
             "Velmeshev_Adults" = "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/Seurat_UCSC/",
             "DISCO" = "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/")

ds_names <- c("Velmeshev_3rd_trimester"= "Velmeshev_2022_3rd_trimester",
              "Velmeshev_0_1_years"= "Velmeshev_2022_0_1_years",
              "Velmeshev_1_2_years"= "Velmeshev_2022_1_2_years",
              "Velmeshev_2_4_years"= "Velmeshev_2022_2_4_years",
              "Velmeshev_Adults"= "Velmeshev_2022_Adult",
              "DISCO" = "brainV1.0_all_FM_filt")


xist_ls <- list()
cell_info <- list()
for (id in names(ds_paths)) {
  print(id)
  if (id == "DISCO") {
    id_seurat <- readRDS(paste0(ds_paths[[id]], ds_names[[id]], ".rds"))
    xist_ls <- append(xist_ls, list(data.frame("XIST"=GetAssayData(id_seurat[["RNA"]], slot="data")["XIST",])))
    cell_info <- append(cell_info, list(read.csv(paste0(ds_paths[[id]][1], "DEGs_common/cell_info.csv")))) 
    rm(id_seurat)
  } else {
    id_seurat <- readRDS(paste0(ds_paths[[id]], ds_names[[id]], ".rds"))
    xist_ls <- append(xist_ls, list(data.frame("XIST"=GetAssayData(id_seurat[["RNA"]], slot="data")["XIST",])))
    cell_info <- append(cell_info, list(read.csv(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/", ds_names[[id]], "/cell_info_", ds_names[[id]], ".csv"))))
    rm(id_seurat)
  }
}
names(xist_ls) <- names(ds_paths)
names(cell_info) <- names(ds_paths)

# 2. Manually add Velmeshev_10_20 years
xist_ls <- append(xist_ls, list("Velmeshev_10_20_years"=read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_2022_10_20_years_XIST.csv", col.names = "XIST")))
cell_info <- append(cell_info, list("Velmeshev_10_20_years"=read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_10_20_years/cell_info_Velmeshev_2022_10_20_years.csv")))

# 3. Merge XIST in df and formats the information
xist_df <- do.call(rbind, xist_ls)
xist_df <- cbind("id"=rownames(xist_df), xist_df)
rownames(xist_df) <- NULL
ids_split <-regmatches(xist_df$id, regexpr("\\.", xist_df$id), invert = TRUE)
ids_split <- as.data.frame(do.call(rbind, ids_split))
xist_df$group <- ids_split$V1
xist_df$sample_id <- ids_split$V2
xist_df <- cbind(xist_df[, c(3:4)], xist_df[, c(2)])

names(xist_df)[names(xist_df)=="xist_df[, c(2)]"] <- "XIST"

# 4. Adds the metadata about the cells
xist_df$info <- rep(NA, nrow(xist_df))
for (group in unique(xist_df$group)) {
  print(group)
  for (id in unique(xist_df[which(xist_df$group==group), "sample_id"])) {
    xist_df[which(xist_df$group==group & xist_df$sample_id==id), "info"] <- cell_info[[group]][which(cell_info[[group]]$cell_id==id), "og_group"]
  }
}

# 5. Uniforms the sex metadata
xist_df$info <- str_replace_all(xist_df$info, c("Female"="F", "Male"="M"))

# 6. Formats the ids so they can be split into barcodes and samples
xist_df$new_ids <- xist_df$sample_id
xist_df$new_ids <- str_replace_all(xist_df$new_ids, c("-1_"="-1/", "-1--"="-1/"))
xist_df <- separate(xist_df, new_ids, into=c("barcodes", "samples"), sep="/")

# 7. Creates sex column
xist_df$sex <- rep(NA, nrow(xist_df))
xist_df[grep("M_", xist_df$info), "sex"] <- "M"
xist_df[grep("F_", xist_df$info), "sex"] <- "F"

# 8. XIST RNA expression plot
groups_order <- c("Velmeshev_3rd_trimester", "Velmeshev_0_1_years",     "Velmeshev_1_2_years",     "Velmeshev_2_4_years",    
                  "Velmeshev_10_20_years" , "Velmeshev_Adults",        "DISCO")
xist_df$group <- factor(xist_df$group, groups_order)
xist_df$samples <- factor(xist_df$samples)
xist_df <- xist_df[order(xist_df$samples), ]

ggplot(xist_df, aes(factor(samples, levels = levels(xist_df$samples)), XIST, fill=sex)) +
  geom_violin() +
  facet_grid(group ~ sex, scales = "free") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"))


# 9. Number of samples
cell_info <- append(cell_info, list("Velmeshev_2nd_trimester"=read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2nd_trimester/cell_info_Velmeshev_2022_2nd_trimester.csv")))
