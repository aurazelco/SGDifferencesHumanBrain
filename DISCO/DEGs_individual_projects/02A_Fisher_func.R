# 0. Import libraries
library(reshape)
library(ggplot2)

# 1. Import data
ImportSignDE <- function(main_dir, ext, row_col) {
  if (missing(ext)) {
    deg_files <- list.files(path = main_dir, pattern = "chr.csv$",full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  else {
    deg_files <- list.files(path = main_dir, pattern = paste0("chr.csv",ext,"$"),full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  names_deg <- list.files(path = main_dir, pattern = "chr.csv$",full.names = FALSE)
  names(deg) <- substr(names_deg, 1, nchar(names_deg)-4)
  return(deg)
}

# 2. Count F and M intersected genes
ReadRawData <- function(main_dir, dis_type, sex, ext, row_col) {
  path <- paste0(main_dir, dis_type, "/01A_DEGs")
  #path <- paste0(main_dir, dis_type, "/01B_num_DEGs")
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
  return(df_sex)
}

# 2. Count F and M intersected genes
ReadRawData2 <- function(main_dir, dis_type, sex, ext, row_col) {
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
  return(df_sex)
}

# 3. Count Unique DEGs 
CountGenes <- function(df_sex) {
  unique_degs <- vector()
  for (i in names(df_sex)) {
    unique_degs <- c(unique_degs, rownames(df_sex[[i]]))
  }
  return(length(unique(unique_degs)))
}

CountGenes2 <- function(df_sex) {
  unique_degs <- vector()
  for (i in names(df_sex)) {
    unique_degs <- c(unique_degs, rownames(df_sex[[i]]))
  }
  return(length(unique(unique_degs)))
}

# 4. Calculate the total number of DEGs
CalcTotGenes <- function(main_dir, dis_type, sex, ext, row_col) {
  df_sex <- ReadRawData(main_dir, dis_type, sex, ext, row_col)
  uniq_sex <- CountGenes(df_sex)
  return(uniq_sex)
}

CalcTotGenes2 <- function(main_dir, dis_type, sex, ext, row_col) {
  df_sex <- ReadRawData2(main_dir, dis_type, sex, ext, row_col)
  uniq_sex <- CountGenes(df_sex)
  return(uniq_sex)
}

# 5. Retrieve num of X and Y-DEGs
XYdeg <- function(main_dir, dis_type, sex, ext, row_col) {
  path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  XY <- ImportSignDE(path, ext, row_col)
  col_factors <- c("ct", "chr")
  sex_index <- paste0(sex, "_num_chr")
  XY[[sex_index]][col_factors] <- lapply(XY[[sex_index]][col_factors], as.factor) 
  df_sex <- data.frame(levels(XY[[sex_index]]$ct))
  df_sex_col <- vector()
  for (chr in levels(XY[[sex_index]]$chr)) {
    df_sex_col <- c(df_sex_col, chr)
    df_sex <- cbind(df_sex, XY[[sex_index]][which(XY[[sex_index]]$chr==chr), "count"])
  }
  rownames(df_sex) <- df_sex[,1]
  df_sex[,1] <- NULL
  colnames(df_sex) <- df_sex_col
  return(df_sex)
}

# 6. Calculate Hypergeometric Distribution
CalcHyperGeom <- function(main_dir, dis_type, tot_genes, df_sex, chr_genes, chr, ext, row_col) {
  pval_res <- vector()
  for (ct in rownames(df_sex)) {
    pval <- phyper(
      df_sex[ct, chr] - 1,
      rowSums(df_sex[ct, ]),
      tot_genes - rowSums(df_sex[ct, ]),
      chr_genes,
      lower.tail= FALSE
    )
    pval_res <- c(pval_res, pval)
  }
  #result_df <- data.frame(rownames(df_sex), pval_res)
  #colnames(result_df) <- c("ct", "pval_hypergeom")
  #result_df$pval_sign <- result_df$pval_hypergeom<0.05
  #if (any(result_df$pval_sign) == TRUE) {
  #  result_df <- result_df[which(result_df$pval_sign==TRUE),]
  #}
  #return(result_df)
  return(pval_res)
}

CalcHyperGeom2 <- function(main_dir, dis_type, tot_genes, df_sex, chr_genes, chr, ext, row_col) {
  pval_res <- vector()
  for (ct in rownames(df_sex)) {
    pval <- phyper(
      df_sex[ct, chr] - 1,
      rowSums(df_sex[ct, ]),
      tot_genes - chr_genes,
      chr_genes,
      lower.tail= FALSE
    )
    pval_res <- c(pval_res, pval)
  }
  return(pval_res)
}

# https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
# depletion
# phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE)
# enrichment -> using this one
# phyper(Overlap-1, group2, Total-group2, group1,lower.tail= FALSE)
# group1 : X_chr_genes
# overlap: df_sex[i, "X"]
# group2: rowSums(df_sex[i,])
# Total: tot sign genes (calculated by CalcTotGenes)

# 7. Calculate HyperGeom for both sexes and chromosomes genes
HyperGeomDEGs <- function(main_dir, dis_type, sex, tot_genes, X_chr, Y_chr, ext, row_col) {
  out_path <- paste0(main_dir, dis_type, "/outputs/02A_Fisher_sex_genes/")
  dir.create(out_path, showWarnings = FALSE, recursive = T)
  df_sex <- XYdeg(main_dir, dis_type, sex, ext, row_col)
  df_sex$X_enriched_pval <- CalcHyperGeom(main_dir, dis_type, tot_genes, df_sex, X_chr, "X", ext, row_col)
  df_sex$Y_enriched_pval <- CalcHyperGeom(main_dir, dis_type, tot_genes, df_sex, Y_chr, "Y", ext, row_col)
  write.csv(df_sex, paste0(out_path, sex, "_Fisher_results.csv"))
  return(df_sex)
}

HyperGeomDEGs2 <- function(main_dir, dis_type, sex, tot_genes, X_chr, Y_chr, ext, row_col) {
  out_path <- paste0(main_dir, dis_type, "/outputs/02A_Fisher_sex_genes/")
  dir.create(out_path, showWarnings = FALSE, recursive = T)
  df_sex <- XYdeg(main_dir, dis_type, sex, ext, row_col)
  df_sex$X_enriched_pval <- CalcHyperGeom2(main_dir, dis_type, tot_genes, df_sex, X_chr, "X", ext, row_col)
  df_sex$Y_enriched_pval <- CalcHyperGeom2(main_dir, dis_type, tot_genes, df_sex, Y_chr, "Y", ext, row_col)
  write.csv(df_sex, paste0(out_path, sex, "_Fisher_results.csv"))
  return(df_sex)
}

# 8. Calculate all and plot results
PlotHyperGeom <- function(main_dir, dis_type, df_sex, sex) {
  out_path <- paste0(main_dir, dis_type, "/outputs/02A_Fisher_sex_genes/")
  dir.create(out_path, showWarnings = FALSE, recursive = T)
  df_melt <- melt.data.frame(df_sex[, c(4:5)])   
  df_melt$ct <- as.factor(rep(rownames(df_sex), 2))
  colnames(df_melt) <- c("chr", "pval", "ct")
  df_melt$chr <- as.factor(df_melt$chr)
  if (any(df_melt$pval< 0.05)) {
    pdf(paste0(out_path, sex, "_pval.pdf"))
    print(ggplot(df_melt[which(df_melt$pval < 0.05), ], aes(ct, chr)) +
            geom_tile(aes(fill=pval)) +
            labs(x="Cell types", y="Chromosome-enriched genes", fill="P-values", main = dis_type) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  axis.title.x = element_blank(),
                  axis.text.x = element_text(size=8, colour = "black", angle = 90, vjust = 0.7, hjust=0.5),
                  axis.ticks.x=element_blank(),
                  axis.title.y = element_text(size=12, face="bold", colour = "black"),
                  legend.position = "bottom", 
                  legend.title = element_text(size=12, face="bold", colour = "black"))
    )
    dev.off()
  }
}

PlotHyperGeom2 <- function(main_dir, dis_type, df_sex, sex) {
  out_path <- paste0(main_dir, dis_type, "/outputs/02A_Fisher_sex_genes/")
  dir.create(out_path, showWarnings = FALSE, recursive = T)
  df_melt <- melt.data.frame(df_sex[, c(4:5)])   
  df_melt$ct <- as.factor(rep(rownames(df_sex), 2))
  colnames(df_melt) <- c("chr", "pval", "ct")
  df_melt$chr <- as.factor(df_melt$chr)
  if (any(df_melt$pval< 0.05)) {
    pdf(paste0(out_path, sex, "_pval_v2.pdf"))
    print(ggplot(df_melt[which(df_melt$pval < 0.05), ], aes(ct, chr)) +
    geom_tile(aes(fill=pval)) +
    labs(x="Cell types", y="Chromosome-enriched genes", fill="P-values", main = dis_type) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=8, colour = "black", angle = 90, vjust = 0.7, hjust=0.5),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
    )
    dev.off()
  }
}
  
# MAIN
SexChr <- function(main_dir, dis_type, tot_genes, X_chr, Y_chr, ext, row_col) {
  df_F <- HyperGeomDEGs(main_dir, dis_type, "F", tot_genes, X_chr, Y_chr, ext, row_col)
  df_M <- HyperGeomDEGs(main_dir, dis_type, "M", tot_genes, X_chr, Y_chr, ext, row_col)
  print("Only significant p-values will be plotted and saved as heatmap")
  PlotHyperGeom(main_dir, dis_type, df_F, "F")
  PlotHyperGeom(main_dir, dis_type, df_M, "M")
}

SexChr2 <- function(main_dir, dis_type, tot_genes, X_chr, Y_chr, ext, row_col) {
  df_F <- HyperGeomDEGs2(main_dir, dis_type, "F", tot_genes, X_chr, Y_chr, ext, row_col)
  df_M <- HyperGeomDEGs2(main_dir, dis_type, "M", tot_genes, X_chr, Y_chr, ext, row_col)
  print("Only significant p-values will be plotted and saved as heatmap")
  PlotHyperGeom2(main_dir, dis_type, df_F, "F")
  PlotHyperGeom2(main_dir, dis_type, df_M, "M")
}
