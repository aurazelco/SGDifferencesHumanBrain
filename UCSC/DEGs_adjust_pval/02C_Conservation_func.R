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
ReadRawData <- function(main_dir, sex, ext, row_col) {
  path <- paste0(main_dir,  "/01B_num_DEGs")
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
  names(df_sex) <- sapply(1:length(names(df_sex)), function(i) str_replace(names(df_sex)[i], paste0("_", sex, "_filt"), ""))
  return(df_sex)
}


# 3. Retrieve all genes expressed in the cell types in a specific disease - from DISCO
AllGenes <- function(all_df, ct_ordered) {
  all_df <- subset(all_df, subset = ct %in% ct_ordered)
  all_df$ct <- droplevels(all_df$ct)
  return(all_df)
}

# 4. Calculate Frction of Sex-biased conserved genes
SexFrac <- function(cons_filt, df_sex, sex) {
  sex_cons <- sapply(1:length(names(df_sex)), function(x) length(intersect(rownames(df_sex[[x]]), cons_filt$gene_name)))
  fr_sex <- sapply(1:length(names(df_sex)), function(x) sex_cons[x] / nrow(df_sex[[x]]))
  return(fr_sex)
}

# 5. Calculate Fraction of all DEGs
AllFrac <- function(cons_filt, genes_df, ct_names, sex) {
  fr_all <- vector()
  for (ct in ct_names) {
    ct_genes <- genes_df[which(genes_df$ct==ct & genes_df$sex==sex), "genes"]
    fr_ct <- length(intersect(ct_genes, cons_filt$gene_name)) / length(ct_genes)
    fr_all <- c(fr_all, fr_ct)
  }
  return(fr_all)
}

# 6. Create new Df
DfFrac <- function(main_dir,  cons_df, threshold, out_name, all_df, ct_ordered) {
  df_F <- ReadRawData(main_dir,  "F")
  df_M <- ReadRawData(main_dir,  "M")
  genes_df <- AllGenes(all_df, ct_ordered)
  if (out_name == "Primates") {
    cons_filt <- subset(cons_df, rowSums(cons_df[, c(5:10)])>=threshold)
  } else if (out_name == "SAGD") {
    cons_filt <- subset(cons_df, rowSums(cons_df[, c(2:22)])>=threshold)
  } else if (out_name == "ENSEMBL") {
    cons_filt <- subset(cons_df, rowSums(cons_df[, c(2:ncol(cons_df))])>=threshold)
  } 
  fr_all_F <- AllFrac(cons_filt, genes_df, names(df_F), "F")
  fr_all_M <- AllFrac(cons_filt, genes_df, names(df_M), "M")
  fr_F <- SexFrac(cons_filt, df_F, "F")
  fr_M <- SexFrac(cons_filt, df_M, "M")
  df_frac <- data.frame(c(rep("F", length(names(df_F))), rep("M", length(names(df_M)))),
                          c(names(df_F), names(df_M)), 
                          c(fr_all_F, fr_all_M),
                          c(fr_F, fr_M)
                          )
  colnames(df_frac) <- c("sex", "ct", "All", "DEG_fraction")
  df_frac$ct <- as.factor(df_frac$ct)
  dir.create(paste(main_dir,  "02C_Conservation", sep="/"), showWarnings = FALSE)
  write.csv(df_frac, paste0(main_dir, "/02C_Conservation/", out_name, "_fraction_in_", threshold,  "_species.csv"), row.names = F)
  df_frac <- melt(df_frac)
  names(df_frac)[names(df_frac) == 'variable'] <- 'group'
  names(df_frac)[names(df_frac) == 'value'] <- 'fractions'
  df_frac$group <- paste(df_frac$sex, df_frac$group, sep="_")
  return(df_frac)
}

# 7. Plot
PlotFrac <- function(main_dir,  df_frac, threshold, out_name, ct_ordered) {
  dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(df_frac$ct))]
  df_frac$ct <- factor(df_frac$ct, dis_ct_ordered)
  df_frac <- df_frac[order(df_frac$ct), ]
  pdf(paste0(main_dir,  "/02C_Conservation/", out_name, "_fraction_in_", threshold,  "_species.pdf"))  
  print(
  ggplot(df_frac, aes(ct, fractions, fill=group)) +
    geom_bar(stat='identity', position='dodge', color="black") + 
    labs(title=paste0(out_name, " - conserved in at least ",threshold, " species"), x="Cell types", y=paste0("Fraction of conserved genes"), fill="Groups") +
    scale_fill_discrete(labels=c("All Female Genes", "Female-biased genes", "All Male Genes", "Male-biased genes")) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  )
  dev.off()
}

# 8. MAIN
ConservedFractions <- function(main_dir,  cons_df, threshold, out_name, all_df, ct_ordered) {
  df_frac <- DfFrac(main_dir,  cons_df, threshold, out_name, all_df, ct_ordered)
  PlotFrac(main_dir,  df_frac, threshold, out_name, ct_ordered)
}