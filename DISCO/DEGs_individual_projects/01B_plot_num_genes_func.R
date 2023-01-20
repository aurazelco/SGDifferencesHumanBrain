################ USER-DEFINED FUNCTIONS

# 0. Import libraries
library(reshape2)
library(ggplot2)
library(stringr)

# 1. Import data
ImportDE <- function(path, ext, row_col) {
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

# 2. Function to filter out ns genes and too low FC, and order based on FC 
Filter_gene <- function( order.gene.df, pval,
                         FC) {
  logFC <- log2(1/FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  
  #If there are sig genes add index number for each sig gene
  
  if(nrow(gene.sig) > 0) {
    gene.sig$index <- seq.int(nrow(gene.sig))
  }
  #print(nrow(gene.sig))
  return(gene.sig)
}

# 3. Count DEGs before intersection and plot
CreateDEGdf <- function(deg_list) {
  proj <- vector()
  sex <- vector()
  ct  <- vector()
  deg_count <- vector()
  for (ct_name in names(deg_list)) {
    ct <- c(ct, rep(ct_name, 2))
    sex <- c(sex, c("F", "M"))
    for (sex_df in names(deg_list[[ct_name]])) {
      deg_count <- c(deg_count, nrow(deg_list[[ct_name]][[sex_df]]))
    }
  }
  df_count <- data.frame(ct, sex, deg_count)
  col_factors <- c("ct", "sex")
  df_count[col_factors] <- lapply(df_count[col_factors], as.factor)  
  return(df_count)
}

# 4. Plot DEGs before intersection
PlotDEGs <- function(main_dir, df_count, ct_ordered) {
  dis_ct_order <- ct_ordered[which(ct_ordered %in% levels(df_count$ct))]
  df_count$ct <- factor(df_count$ct, dis_ct_order)
  df_count <- df_count[order(df_count$ct), ]
  pdf(paste0(main_dir, "num_DEGs.pdf"))
  print(
    ggplot(df_count, aes(ct, deg_count, fill=sex)) +
      geom_bar(stat="identity", position = "dodge") +
      labs(x="Cell types", y="Number of DEGs", fill="Sex") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, angle = 90, colour = "black"),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.title = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom")
  )
  dev.off()
}

# 5. Extract gene lists and number of common DEGs
CountDEG <- function(main_dir, dis_type, pval, FC, ct_ordered) {
  #print("The script does not check if a folder has only one project - this may raise errors")
  out_path <- paste0(main_dir, dis_type, "/outputs/01B_num_DEGs/")
  dir.create(out_path, showWarnings = FALSE)
  path <- (paste0(main_dir, dis_type, "/outputs/01A_DEGs/"))
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  celltypes <- rep(sub_ct, each=2)
  sex <- rep(c("F", "M"), length.out = length(sub_ct))
  df <- data.frame(cbind(celltypes, sex))
  num_genes <- vector()
  ct_qc <- vector()
  deg_list <- list()
  for (ct in 1:length(sub_ct)) {
    all_DEGs <- ImportDE(paste0(path, sub_ct[ct]))
    filt_DEGs <- lapply(all_DEGs, function(x) Filter_gene(x, pval, FC))
    dir.create( paste0(out_path, sub_ct[ct]), recursive = T, showWarnings = F)
    lapply(1:length(names(filt_DEGs)), function(x) write.csv(filt_DEGs[[x]],
                                                             paste0(out_path, sub_ct[ct], "/", names(filt_DEGs)[x], "_filt.csv")
    ))
    deg_list <- append(deg_list, list(filt_DEGs))
  }
  names(deg_list) <- sub_ct
  deg_df <- CreateDEGdf(deg_list)
  PlotDEGs(out_path, deg_df, ct_ordered)
}

