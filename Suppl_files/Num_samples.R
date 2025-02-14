# Author: Aura Zelco
# Brief description:
  # This script is used for calculating the number of samples used per study
  # Brief procedure:
    # 1. it imports the samples information for all datasets
    # 2. it calculates the number of samples per dataset based on the barcodes
    # 3. it adds the study ID
    # 4. sums up the number of samples per study and saves it as a LaTEX table
# Documentation abbreviations:
  # F and M: females and males
  # df: dataframe
  # ds: dataset

#---------------------------------------------------------------------------------------------------

# 0. Import Libraries

library(stringr)
library(tidyr)
library(ggplot2)
library(xtable)

# 1. Imports cell information for all datasets

main <- "data/"

ds_paths <- c("Velmeshev_3rd_trimester"= "UCSC/Seurat_UCSC/",
              "Velmeshev_0_1_years"="UCSC/Seurat_UCSC/",
              "Velmeshev_1_2_years"= "UCSC/Seurat_UCSC/",
              "Velmeshev_2_4_years"= "UCSC/Seurat_UCSC/",
              "Velmeshev_Adults" = "UCSC/Seurat_UCSC/",
              "DISCO" = "DISCOv1.0/")

ds_names <- c("Velmeshev_3rd_trimester"= "Velmeshev_2022_3rd_trimester",
              "Velmeshev_0_1_years"= "Velmeshev_2022_0_1_years",
              "Velmeshev_1_2_years"= "Velmeshev_2022_1_2_years",
              "Velmeshev_2_4_years"= "Velmeshev_2022_2_4_years",
              "Velmeshev_Adults"= "Velmeshev_2022_Adult",
              "DISCO" = "brainV1.0_all_FM_filt")


cell_info <- list()
for (id in names(ds_paths)) {
  print(id)
  if (id == "DISCO") {
    cell_info <- append(cell_info, list(read.csv(paste0(main, ds_paths[[id]][1], "DEGs_common/cell_info.csv")))) 
  } else {
    cell_info <- append(cell_info, list(read.csv(paste0(main, "UCSC/DEGs_adjust_pval/", ds_names[[id]], "/cell_info_", ds_names[[id]], ".csv"))))
  }
}
names(cell_info) <- names(ds_paths)

cell_info <- append(cell_info, list("Velmeshev_2nd_trimester"=read.csv(paste0(main, "UCSC/DEGs_2nd_trimester/Velmeshev_2022_2nd_trimester/cell_info_Velmeshev_2022_2nd_trimester.csv"))))

# 2. Creates a dataframe with the number of samples in each dataset based on the barcode

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


num_samples_simplified$study <- paste(str_split(num_samples_simplified$new_group, "_", simplify=T)[,1], num_samples_simplified$sex, sep="--")

# 3. Sums the number of samples in each study

study_samples <- c()
for (study_id in unique(num_samples_simplified$study)) {
  study_samples <- c(study_samples, sum(num_samples_simplified[which(num_samples_simplified$study==study_id), "samples_count"]))
}
study_df <- data.frame(unique(num_samples_simplified$study), study_samples)
colnames(study_df) <- c("study_id", "num_samples")
study_df <- separate(data = study_df, col = "study_id", into = c("study", "sex"), sep="--")
study_df$sex <- as.factor(study_df$sex)
study_df <- study_df[order(study_df$sex),]
study_df <- study_df[order(study_df$study),]

# 4. Prints the resulting table into a LaTEX table

print(xtable(study_df, type = "latex"), file = "Samples_table.tex")