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
    ct <- c(ct, rep(ct_name, length(deg_list[[ct_name]]$F)*2))
    projs <- sapply(1:length(names(deg_list[[ct_name]]$F)), function(i) str_remove(names(deg_list[[ct_name]]$F)[i], "_F"))
    proj <- c(proj, rep(projs, 2))
    sex <- c(sex, rep(c("F", "M"), each = length(deg_list[[ct_name]]$F)))
    deg_count <- c(deg_count, sapply(1:length(names(deg_list[[ct_name]]$F)), function(i) nrow(deg_list[[ct_name]]$F[[i]])))
    deg_count <- c(deg_count, sapply(1:length(names(deg_list[[ct_name]]$M)), function(i) nrow(deg_list[[ct_name]]$M[[i]])))
  }
  df_count <- data.frame(ct, proj, sex, deg_count)
  col_factors <- c("ct", "proj", "sex")
  df_count[col_factors] <- lapply(df_count[col_factors], as.factor)  
  return(df_count)
}

# 4. Plot DEGs before intersection
PlotDEGs <- function(main, dis, df_count, ct_ordered) {
  dis_ct_order <- ct_ordered[which(ct_ordered %in% levels(df_count$ct))]
  df_count$ct <- factor(df_count$ct, dis_ct_order)
  df_count <- df_count[order(df_count$ct), ]
  pdf(paste0(main, dis,  "/01B_num_DEGs/num_genes_per_ct_before_intersection.pdf"),
      width = 8)
  print(
  ggplot(df_count, aes(ct, deg_count, fill=ct)) +
    geom_bar(stat="identity") +
    facet_wrap(~sex*proj, scales="free") +
    labs(y="Number of DEGs", fill="") +
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

# 5. Extract gene lists and number of common DEGs
IntersectDEG <- function(main, dis, pval, FC, ct_ordered) {
  #print("The script does not check if a folder has only one project - this may raise errors")
  dir.create(paste(main, dis, "01B_num_DEGs", sep="/"), showWarnings = FALSE)
  path <- (paste(main, dis, "01A_DEGs", sep="/"))
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  celltypes <- rep(sub_ct, each=2)
  sex <- rep(c("F", "M"), length.out = length(sub_ct))
  df <- data.frame(cbind(celltypes, sex))
  num_genes <- vector()
  ct_qc <- vector()
  deg_list <- list()
  num_proj <- vector()
  for (ct in 1:length(sub_ct)) {
    all_DEGs <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    filt_DEGs <- lapply(all_DEGs, function(x) Filter_gene(x, pval, FC))
    df_F <- list()
    df_M <- list()
    for (i in names(filt_DEGs)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, all_DEGs[i])
      } else {
        df_M <- append(df_M, all_DEGs[i])
      }
    }
    ct_list <- list("F" = df_F, "M" = df_M)
    deg_list <- append(deg_list, list(ct_list))
    if (length(df_F)>1) {
      dir.create(paste(main, dis, "01B_num_DEGs", sub_ct[ct], sep="/"), showWarnings = FALSE)
      num_proj <- c(num_proj, length(df_F))
      if (length(df_F)==3) {
        int_list <- list(Reduce(intersect, list(rownames(df_F[[1]]), rownames(df_F[[2]]), rownames(df_F[[3]]))),
                         Reduce(intersect, list(rownames(df_M[[1]]), rownames(df_M[[2]]), rownames(df_M[[3]]))))
        names(int_list) <- c("F", "M")
        lapply(1:length(int_list), function(i) write.csv(int_list[i], 
                                                         file = paste0(main, "/", dis, "/01B_num_DEGs/", sub_ct[ct], "/", names(int_list[i]), "_intersected_genes.csv"),
                                                         row.names = TRUE))
      } else {
        int_list <- list(intersect(rownames(df_F[[1]]), rownames(df_F[[2]]))
                         ,intersect(rownames(df_M[[1]]), rownames(df_M[[2]])))
        names(int_list) <- c("F", "M")
        lapply(1:length(int_list), function(i) write.csv(int_list[i], 
                                                         file = paste0(main, "/", dis,  "/01B_num_DEGs/", sub_ct[ct],"/", names(int_list[i]), "_intersected_genes.csv"),
                                                         row.names = TRUE))
      }
      num_genes <- c(num_genes, length(int_list[[1]]))
      num_genes <- c(num_genes, length(int_list[[2]]))
      ct_qc <- c(ct_qc, sub_ct[ct])
      ct_qc <- c(ct_qc, sub_ct[ct])
    } else {
      dir.create(paste(main, dis, "01B_num_DEGs", sub_ct[ct], sep="/"), showWarnings = FALSE)
      num_proj <- 1
      int_list <- list(rownames(df_F[[1]])
                       ,rownames(df_M[[1]]))
      names(int_list) <- c("F", "M")
      lapply(1:length(int_list), function(i) write.csv(int_list[i], 
                                                       file = paste0(main, "/", dis,  "/01B_num_DEGs/", sub_ct[ct],"/", names(int_list[i]), "_intersected_genes.csv"),
                                                       row.names = TRUE))
      num_genes <- c(num_genes, nrow(df_F[[1]]))
      num_genes <- c(num_genes, nrow(df_M[[1]]))
      ct_qc <- c(ct_qc, sub_ct[ct])
      ct_qc <- c(ct_qc, sub_ct[ct])
    }
  }
  names(deg_list) <- sub_ct
  deg_df <- CreateDEGdf(deg_list)
  PlotDEGs(main, dis, deg_df, ct_ordered)
  if (identical(df$celltypes, ct_qc)) {
    df$counts <- num_genes
    col_factors <- c("celltypes", "sex")
    df[col_factors] <- lapply(df[col_factors], as.factor) 
  } else {
    print("The cell types did not match when intersected; please check the input data")
  }
  return(list(df, num_proj))
}

# 6. Plot number of intersected DEGs
PlotIntDEGs <- function(main, dis, df, num_p, ct_ordered) {
  dis_ct_order <- ct_ordered[which(ct_ordered %in% levels(df$celltypes))]
  df$celltypes <- factor(df$celltypes, dis_ct_order)
  df <- df[order(df$celltypes), ]
  pdf(paste0(main, "/", dis,  "/01B_num_DEGs/num_genes_per_ct.pdf"))
  print(ggplot(df, aes(celltypes, counts, fill=sex)) +
          geom_bar(stat="identity", position="dodge") +
          labs(x="", y="Intersected number of DEGs", fill="Sex", main = dis) +
          #annotate("text", x = seq(1,length(levels(df$celltypes))), y = -1, label = num_p) +
          annotate("text", x = seq(1,length(levels(df$celltypes))), y =(df$counts[c(TRUE, FALSE)] + 10) , label = num_p, colour = "blue") +
          annotate("text", x = length(levels(df$celltypes))*0.75, y =max(df$counts) + 30 , label = "Number of intersected projects", colour = "blue") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black")))
  dev.off()
}

# 