# 0. Import libraries
#library(readxl)
library(stringr)
library(reshape)
library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Import data
ImportSignDE <- function(main_dir, ext, row_col) {
  if (missing(ext)) {
    deg_files <- list.files(path = main_dir, pattern = "\\.csv$",full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  else {
    deg_files <- list.files(path = main_dir, pattern = paste0("\\.",ext,"$"),full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  names_deg <- list.files(path = main_dir, pattern = "\\.csv$",full.names = FALSE)
  names(deg) <- substr(names_deg, 1, nchar(names_deg)-4)
  return(deg)
}

# 2. Retrieve DEGs
ReadRawData <- function(main_dir, dis_type, sex, ext, row_col) {
  path <- paste0(main_dir, dis_type, "/01B_num_DEGs")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_sex <- list()
  df_names <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportSignDE(paste(path, sub_ct[ct], sep="/"))
    for (i in names(deg)) {
      if (grepl(sex, i, fixed=TRUE)){
        df_names <- c(df_names, paste(sub_ct[ct], i, sep="_"))
        df_sex <- append(df_sex, list(deg[[i]]))
      }
    }
  }
  names(df_sex) <- df_names
  names(df_sex) <- sapply(1:length(names(df_sex)), function(i) str_replace(names(df_sex)[i], paste0("_", sex, "_intersected_genes"), ""))
  return(df_sex)
}

# 3. Retrieve all DEGs file
AllDEGs <- function(all_degs, dis_type) {
  all_degs[, "gene"] <- NULL
  names(all_degs)[names(all_degs) == 'X'] <- 'gene_name'
  all_degs$cluster <- str_replace_all(all_degs$cluster, "/", "_")
  all_degs$cluster <- as.factor(all_degs$cluster)
  if (dis_type == "Normal") {
    all_degs <- subset(all_degs, subset = cluster !=  c("T"))
  } else if (dis_type == "Alzheimer's disease") {
    ad_remove <- list.dirs("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/Alzheimer's disease/01A_only_1_project", recursive=FALSE, full.names = FALSE)
    `%!in%` <- Negate(`%in%`)
    all_degs <- subset(all_degs, subset = cluster %!in% ad_remove)
  }
  all_degs$cluster <- droplevels(all_degs$cluster)
  return(all_degs)
}

# 4. Calculate Frction of Sex-biased conserved genes
SexFrac <- function(cons_filt, df_sex, sex) {
  sex_cons <- sapply(1:length(names(df_sex)), function(x) length(intersect(df_sex[[x]][[sex]], cons_filt$gene_name)))
  fr_sex <- sapply(1:length(names(df_sex)), function(x) sex_cons[x] / length(df_sex[[x]][[sex]]))
  return(fr_sex)
}

# 5. Calculate Fraction of all DEGs
AllFrac <- function(cons_filt, all_degs, ct_names) {
  fr_all <- vector()
  for (ct in ct_names) {
    ct_deg <- subset(all_degs, subset = cluster == ct)
    ct_deg$cluster <- droplevels(ct_deg$cluster)
    fr_all <- c(fr_all, (length(intersect(ct_deg$gene_name, cons_filt$gene_name)) / nrow(ct_deg)))
  }
  print(fr_all)
  return(fr_all)
}

# 6. Create new Df
DfFrac <- function(main_dir, dis_type, cons_df, threshold, out_name, all_degs) {
  df_F <- ReadRawData(main_dir, dis_type, "F")
  df_M <- ReadRawData(main_dir, dis_type, "M")
  all_degs <- AllDEGs(all_degs, dis_type)
  if (out_name == "Primates") {
    cons_filt <- subset(cons_df, rowSums(cons_df[, c(5:10)])>=threshold)
  } else if (out_name == "SAGD") {
    cons_filt <- subset(cons_df, rowSums(cons_df[, c(2:22)])>=threshold)
  }
  if (length(names(df_F)) == length(names(df_M))) {
    #fr_all <- rep(nrow(cons_filt) / nrow(cons_df), length.out = length(names(df_F)))
    fr_all <- AllFrac(cons_filt, all_degs, names(df_F))
  }
  fr_F <- SexFrac(cons_filt, df_F, "F")
  fr_M <- SexFrac(cons_filt, df_M, "M")
  df_frac <- data.frame(names(df_F), fr_all, fr_F, fr_M)
  colnames(df_frac) <- c("ct", "All", "F", "M")
  df_frac$ct <- as.factor(df_frac$ct)
  dir.create(paste(main_dir, dis_type, "02C_Conservation", sep="/"), showWarnings = FALSE)
  write.csv(df_frac, paste0(main_dir, "/", dis_type, "/02C_Conservation/", out_name, ".csv"))
  df_frac <- melt(df_frac, variable_name = "group")
  names(df_frac)[names(df_frac) == 'value'] <- 'fractions'
  return(df_frac)
}

# 7. Plot
PlotFrac <- function(main_dir, dis_type, df_frac, threshold, out_name) {
  pdf(paste0(main_dir, "/", dis_type, "/02C_Conservation/", out_name, "_fraction_in_", threshold,  "_species.pdf"))  
  print(
  ggplot(df_frac, aes(ct, fractions, fill=group)) +
    geom_bar(stat='identity', position='dodge') + 
    labs(title=paste0(out_name, " - conserved in at least ",threshold, " species"), x="Cell types", y=paste0("Fraction of conserved genes"), fill="Groups") +
    scale_fill_discrete(labels=c("All", "Female-biased genes", "Male-biased genes")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}

# 8. MAIN
ConservedFractions <- function(main_dir, dis_type, cons_df, threshold, out_name, all_degs) {
  df_frac <- DfFrac(main_dir, dis_type, cons_df, threshold, out_name, all_degs)
  PlotFrac(main_dir, dis_type, df_frac, threshold, out_name)
}