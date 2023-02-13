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
    #id_seurat <- readRDS(paste0(ds_paths[[id]], ds_names[[id]], ".rds"))
    #xist_ls <- append(xist_ls, list(data.frame("XIST"=GetAssayData(id_seurat[["RNA"]], slot="data")["XIST",])))
    cell_info <- append(cell_info, list(read.csv(paste0(ds_paths[[id]][1], "DEGs_common/cell_info.csv")))) 
    #rm(id_seurat)
  } else {
    #id_seurat <- readRDS(paste0(ds_paths[[id]], ds_names[[id]], ".rds"))
    #xist_ls <- append(xist_ls, list(data.frame("XIST"=GetAssayData(id_seurat[["RNA"]], slot="data")["XIST",])))
    cell_info <- append(cell_info, list(read.csv(paste0("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/", ds_names[[id]], "/cell_info_", ds_names[[id]], ".csv"))))
    #rm(id_seurat)
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
disco <- as.data.frame(do.call(rbind, str_split(xist_df[which(xist_df$group=="DISCO"), "info"], "_")))

proj_ids <- unique(disco$V1)
dis_ids <- unique(disco$V3)

#for (i in 1:(length(names(cell_info)) -1)) {
#   xist_df[which(xist_df$sample_id %in% cell_info[[i]]$cell_id), "group"] <- names(cell_info)[i]
#}

xist_df$new_group <- xist_df$group
for (proj in proj_ids) {
  for (dis in dis_ids ) {
    for (sex in c("F", "M")) {
      if (any(grep(paste(proj, sex, dis, sep = "_"), xist_df[which(xist_df$group=="DISCO"), "info"]))) {
        print(paste(proj, dis, sep = "_"))
        xist_df[which(grepl(paste(proj, sex, dis, sep = "_"), xist_df$info)), "new_group"] <- paste(proj, dis, sep = "_")
      }
    }
  }
}

groups_order <- c("Velmeshev_3rd_trimester","Velmeshev_0_1_years",          
                  "Velmeshev_1_2_years", "Velmeshev_2_4_years",           
                  "Velmeshev_10_20_years", "Velmeshev_Adults",              
                  "GSE157827_Normal","GSE174367_Normal",              
                  "PRJNA544731_Normal", "GSE157827_Alzheimer's disease",
                  "GSE174367_Alzheimer's disease", "PRJNA544731_Multiple Sclerosis")
xist_df$new_group <- factor(xist_df$new_group, groups_order)
xist_df <- xist_df[order(xist_df$new_group), ]
xist_df$samples <- factor(xist_df$samples, unique(xist_df$samples))
xist_df <- xist_df[order(xist_df$samples), ]


pdf("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Extra_figures/XIST_faceted.pdf", width = 15, height = 20)
print(
  ggplot(xist_df, aes(samples, XIST, fill=sex)) +
    geom_violin() +
    facet_grid(new_group ~ sex, scales = "free") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          strip.text = element_text(size=12, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
)
dev.off()


# 9. Number of samples
cell_info <- append(cell_info, list("Velmeshev_2nd_trimester"=read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/DEGs/Velmeshev_2022_2nd_trimester/cell_info_Velmeshev_2022_2nd_trimester.csv")))
  
num_samples <- cell_info
num_samples$DISCO <- num_samples$DISCO[, c(1,2,5)]
num_samples <- do.call(rbind, num_samples)
num_samples <- cbind("group" = gsub('\\..*', '', rownames(num_samples)), num_samples)
rownames(num_samples) <- NULL
num_samples$X <- NULL

num_samples$new_ids <- num_samples$cell_id
num_samples$new_ids <- str_replace_all(num_samples$new_ids, c("-1_"="-1/", "-1--"="-1/"))
num_samples <- separate(num_samples, new_ids, into=c("barcodes", "samples"), sep="/")
num_samples$sex <- rep(NA, nrow(num_samples))
num_samples$og_group <- str_replace_all(num_samples$og_group, c("Female"="F", "Male"="M"))
num_samples[grep("M_", num_samples$og_group), "sex"] <- "M"
num_samples[grep("F_", num_samples$og_group), "sex"] <- "F"

disco <- as.data.frame(do.call(rbind, str_split(num_samples[which(num_samples$group=="DISCO"), "og_group"], "_")))

proj_ids <- unique(disco$V1)
dis_ids <- unique(disco$V3)


num_samples$new_group <- num_samples$group
for (proj in proj_ids) {
  for (dis in dis_ids ) {
    for (sex in c("F", "M")) {
      if (any(grep(paste(proj, sex, dis, sep = "_"), num_samples[which(num_samples$group=="DISCO"), "og_group"]))) {
        print(paste(proj, dis, sep = "_"))
        num_samples[which(grepl(paste(proj, sex, dis, sep = "_"), num_samples$og_group)), "new_group"] <- paste(proj, dis, sep = "_")
      }
    }
  }
}

num_samples$samples_count <- rep(NA, nrow(num_samples)) 
for (group in unique(num_samples$new_group)) {
  for (sex in c("F", "M")) {
    num_samples[which(num_samples$new_group==group & num_samples$sex==sex), "samples_count"] <- length(unique(num_samples[which(num_samples$new_group==group & num_samples$sex==sex), "samples"]))
  }
}

num_samples_simplified <- num_samples[, c("new_group", "sex", "samples_count")]
num_samples_simplified <- num_samples_simplified[which(!duplicated(num_samples_simplified)),]

sum(num_samples_simplified$samples_count)
sum(num_samples_simplified[which(num_samples_simplified$sex=="F"), "samples_count"])
sum(num_samples_simplified[which(num_samples_simplified$sex=="M"), "samples_count"])


velm_4_10 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_cells_per_age.csv")
velm_4_10 <- subset(velm_4_10, age=="4-10 years")
num_samples_simplified <- rbind(num_samples_simplified,
                                c("Velmeshev_4_10_years", "F", nrow(velm_4_10[which(velm_4_10$sex=="Female"),])),
                                c("Velmeshev_4_10_years", "M", nrow(velm_4_10[which(velm_4_10$sex=="Male"),])))
num_samples_simplified$samples_count <- as.numeric(num_samples_simplified$samples_count)

groups_order <- c("Velmeshev_2nd_trimester", "Velmeshev_3rd_trimester","Velmeshev_0_1_years",          
                  "Velmeshev_1_2_years", "Velmeshev_2_4_years",  "Velmeshev_4_10_years",         
                  "Velmeshev_10_20_years", "Velmeshev_Adults",              
                  "GSE157827_Normal","GSE174367_Normal",              
                  "PRJNA544731_Normal", "GSE157827_Alzheimer's disease",
                  "GSE174367_Alzheimer's disease", "PRJNA544731_Multiple Sclerosis")
num_samples_simplified$new_group <- factor(num_samples_simplified$new_group, groups_order)
num_samples_simplified <- num_samples_simplified[order(num_samples_simplified$new_group), ]

pdf("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Extra_figures/num_samples.pdf")
print(
  ggplot(num_samples_simplified, aes(new_group, samples_count, fill=sex)) +
        geom_bar(stat = "identity", color="black", position = "dodge") +
        geom_hline(yintercept = 3, linetype="dashed") +
        labs(x="Groups", y="Number of samples", fill="Sex") +
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
  )
dev.off()
