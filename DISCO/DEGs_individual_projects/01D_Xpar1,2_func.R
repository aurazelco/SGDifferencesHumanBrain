################ USER-DEFINED FUNCTIONS

# 0. Import libraries
library(ggplot2)
library(dplyr)

# 1. Import data
ImportSignDE <- function(path, ext, row_col) {
  if (missing(ext)) {
    deg_files <- list.files(path = path, pattern = "\\.csv$",full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  else {
    deg_files <- list.files(path = path, pattern = paste0("\\.",ext,"$"),full.names = TRUE)
    if (missing(row_col)) {
      deg <- lapply(deg_files, read.csv, row.names=1)
    }
    else {
      deg <- lapply(deg_files, read.csv, row.names=row_col)
    }
  }
  names_deg <- list.files(path = path, pattern = "\\.csv$",full.names = FALSE)
  names(deg) <- substr(names_deg, 1, nchar(names_deg)-4)
  return(deg)
}

# 2. Intersect with Xpar lists
CompareXpar <- function(gene_list, xpar1, xpar2, sex) {
  xpar1_genes <- intersect(rownames(gene_list), xpar1)
  xpar2_genes <- intersect(rownames(gene_list), xpar2)
  xpar <- list(xpar1_genes, xpar2_genes)
  names(xpar) <- c("Xpar1", "Xpar2")
  xpar <- xpar[lengths(xpar) > 0L]
  return(xpar)
}

# 3. Calculate all Xpars
CalcXpar <- function(genes_df, xpar1, xpar2, sex) {
  Xpar_l <- list()
  Xpar_names <- vector()
  for (i in names(genes_df)) {
    xpar_ct <- CompareXpar(genes_df[[i]], Xpar1_list, Xpar2_list, sex)
    if (length(xpar_ct)>0) {
      Xpar_l <- append(Xpar_l, list(xpar_ct))
      Xpar_names <- c(Xpar_names, i)
    }
  }
  names(Xpar_l) <- Xpar_names
  return(Xpar_l)
}


# 4. Save found Xpar genes as df
SaveXpar <- function(main_dir, dis_type, Xpar_l) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01D_Xpars/")
  dir.create(out_path, showWarnings = F, recursive = T)
  xpar_genes <- vector()
  df_ct <- vector()
  which_xpar <- vector()
  df_sex <- vector()
  for (sex_index in 1:length(Xpar_l)) {
    Xpar_l[[sex_index]] <- Xpar_l[[sex_index]][lengths(Xpar_l[[sex_index]])>0L]
    for (ct_index in 1:length(Xpar_l[[sex_index]])) {
      for (xpar_index in 1:length(Xpar_l[[sex_index]][[ct_index]])) {
        xpar_genes <- c(xpar_genes, Xpar_l[[sex_index]][[ct_index]][[xpar_index]])
        df_ct <- c(df_ct, 
                   rep(names(Xpar_l[[sex_index]][ct_index]), 
                       length.out=length(Xpar_l[[sex_index]][[ct_index]][[xpar_index]])))
        which_xpar <- c(which_xpar, 
                        rep(names(Xpar_l[[sex_index]][[ct_index]][xpar_index]),
                            length.out=length(Xpar_l[[sex_index]][[ct_index]][[xpar_index]])))
        df_sex <- c(df_sex,
                    rep(names(Xpar_l[sex_index]),
                        length.out=length(Xpar_l[[sex_index]][[ct_index]][[xpar_index]])))
      }
    }
  }
  df_xpar <- data.frame(df_sex, df_ct, which_xpar, xpar_genes)
  col_factors <- c("df_sex", "df_ct", "which_xpar")
  df_xpar[col_factors] <- lapply(df_xpar[col_factors], as.factor) 
  write.csv(df_xpar, paste0(out_path, "Xpars_genes.csv"))
  return(df_xpar)
}
  
# 5. Plot
PlotXpar <- function(main_dir, dis_type, df_xpar, ct_ordered) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01D_Xpars/")
  dir.create(out_path, showWarnings = F, recursive = T)
  plt_xpar <- df_xpar[, -c(4)]
  plt_xpar <- as.data.frame(table(plt_xpar))
  names(plt_xpar)[names(plt_xpar) == 'Freq'] <- "num_xpar"
  col_factors <- c("df_sex", "df_ct", "which_xpar")
  plt_xpar[col_factors] <- lapply(plt_xpar[col_factors], as.factor) 
  dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(plt_xpar$df_ct))]
  plt_xpar$df_ct <- factor(plt_xpar$df_ct, dis_ct_ordered)
  plt_xpar <- plt_xpar[order(plt_xpar$df_ct), ]
  pdf(paste0(out_path, "num_xpar.pdf"))
  print(ggplot(plt_xpar, aes(df_ct, num_xpar, fill=df_sex)) +
    geom_bar(stat="identity", position="dodge") +
    facet_wrap(~which_xpar, scales = "free") +
    labs(title=dis_type, x="Cell types", y="Number of Xpar genes", fill="Sex") +
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
  write.csv(plt_xpar, paste0(out_path, "Xpars_genes_table.csv"))
}

# 6. Analyze ct
XparCt <- function(main_dir, dis_type, xpar1, xpar2, ct_ordered, ext, row_col) {
  path <- paste0(main_dir, dis_type, "/outputs/01B_num_DEGs/")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportSignDE(paste0(path, sub_ct[ct]))
    for (i in names(deg)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(deg[[i]]))
        names_F <- c(names_F, sub_ct[ct])
      } else {
        df_M <- append(df_M, list(deg[[i]]))
        names_M <- c(names_M, sub_ct[ct])
      }
    }
  }
  names(df_F) <- names_F
  names(df_M) <- names_M
  XparF <- CalcXpar(df_F, xpar1, xpar2, "F")
  XparM <- CalcXpar(df_M, xpar1, xpar2, "M")
  XparAll <- list(XparF, XparM)
  names(XparAll) <- c("F", "M")
  for (sex in names(XparAll)) {
    if (length(XparAll[[sex]]) == 0) {
      print(paste0("No genes belonging to Xpar1 or Xpar2 in the ", sex, " whole DEG list"))
    }
  }
  if (all(lengths(XparAll)!=0L)) {
    XparAll <- XparAll[lengths(XparAll)>0L]
    XparAll_df <- SaveXpar(main_dir, dis_type, XparAll)
    Xparplot <- PlotXpar(main_dir, dis_type, XparAll_df, ct_ordered)
  } else {
    print("No genes on Xpar1 or Xpar2 in both females and males, all cell types")
  }
}

