################ USER-DEFINED FUNCTIONS

# 0. Import libraries
library(reshape2)
library(ggplot2)
library(stringr)
library(scales)
library(biomaRt)
require(biomaRt)

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

# 2. Function to get chromosome number from gene symbol - Pattama
Annot.chr.name <- function(gene.list){
  
  # define biomart object
  mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",mirror = "useast")
  Annot_idf <- getBM(attributes = c("hgnc_symbol",
                                    "chromosome_name",
                                    "description"),
                     filters = c("hgnc_symbol") ,
                     values=list(gene.list),
                     mart = mart)
  #delete chromosome name with CHR label
  library(stringr)
  Annot_df <- Annot_idf[!str_detect(Annot_idf$chromosome_name,  "CHR"),]
  return(Annot_df)
}

# 3. Map genes from intersected genes against chromosome - Pattama
map_chr <- function(RRA_df,Annot_df){
  map_chr_df <- merge(RRA_df,Annot_df, by.x= "Gene", by.y= "hgnc_symbol")
  return(map_chr_df)
}

# 4. Count the number of genes of X- and Y-chr - Pattama
num_chr2 <- function(map_chr_df){
  chrX <- lapply(map_chr_df, function(x) nrow(x[x$chromosome_name == "X",]))
  chrY <- lapply(map_chr_df, function(x) nrow(x[x$chromosome_name == "Y",]))
  chrA <- lapply(map_chr_df, function(x) nrow(x[ (!x$chromosome_name == "X"|!x$chromosome_name == "Y"),]))
  num_df <- data.frame(X = unlist(chrX), Y = unlist(chrY), Autosome = unlist(chrA))
  return(num_df)
}

# 5. Counts chr and organize df
num_chr_order <- function(df, path, sex) {
  num_df <- num_chr2(df)
  num_df$ct <- as.factor(rownames(num_df))
  num_df <- melt(num_df, id.vars = "ct")
  colnames(num_df) <- c("ct", "chr", "count")
  write.csv(num_df, paste0(path,"/", sex, "_num_chr.csv"))
  return(num_df)
}

# 6. Extract Gene names - X for F-DEGs and Y for for M-DEGs
ExtractSexGenes <- function(chr_sex, chr) {
  sexchr <- lapply(chr_sex, function(x) x[x$chromosome_name == chr,"Gene"])
  genes_names <- unique(unlist(sexchr, use.names = FALSE))
  chr_mtx <- matrix(nrow = length(names(sexchr)), ncol=length(genes_names))
  rownames(chr_mtx) <- names(sexchr)
  colnames(chr_mtx) <- genes_names
  for (ct in rownames(chr_mtx)) {
    for (gene in colnames(chr_mtx)) {
      if (gene %in% sexchr[[ct]]) {
        chr_mtx[ct, gene] <- "y"
      } else {
        chr_mtx[ct, gene] <- "n"
      }
    }
  }
  chr_df <- reshape::melt.matrix(chr_mtx)
  colnames(chr_df) <- c("ct", "gene", "DEG")
  col_factors <- c("ct")
  chr_df[col_factors] <- lapply(chr_df[col_factors], as.factor) 
  return(chr_df)
}

# 7. Extract Gene names for all DEGs
ExtractGenes <- function(chr_list) {
  new_chr_list <- list()
  for (sex in names(chr_list)) {
    df_sex <- do.call(rbind, chr_list[[sex]])
    df_sex$ct <- rownames(df_sex)
    df_sex$ct <- str_remove_all(df_sex$ct, "\\.\\d+")
    rownames(df_sex) <- NULL
    df_sex <- df_sex[, -c(3)]
    df_sex$chr_simplified <- df_sex$chromosome_name
    df_sex$chr_simplified <- str_replace_all(df_sex$chr_simplified, "\\d+", "Autosome")
    df_sex$chr_simplified <- str_replace_all(df_sex$chr_simplified, "MT", "Autosome")
    #col_factors <- c("ct", "Gene", "chromosome_name")
    #df_sex[col_factors] <- lapply(df_sex[col_factors], as.factor) 
    mapped_genes <- unique(df_sex[, c("Gene", "chromosome_name", "chr_simplified")])
    genes_names <- unique(df_sex$Gene)
    chr_mtx <- matrix(nrow = length(unique(df_sex$ct)), ncol=length(genes_names))
    rownames(chr_mtx) <- unique(df_sex$ct)
    colnames(chr_mtx) <- genes_names
    for (i in rownames(chr_mtx)) {
      df_ct <- subset(df_sex, subset = ct == i)
      for (gene in colnames(chr_mtx)) {
        chr_mtx[i, gene] <- ifelse(gene %in% df_ct$Gene, "y", "n")
        #if (gene %in% df_ct$Gene) {
        #  chr_mtx[i, gene] <- "y"
        #} else {
        #  chr_mtx[i, gene] <- "n"
        #}
      }
    }
    #chr_mtx[rev(order(rowSums(chr_mtx))), ]
    chr_df <- reshape::melt.matrix(chr_mtx)
    colnames(chr_df) <- c("ct", "gene", "DEG")
    chr_df$chr_name <- rep(NA, length=nrow(chr_df))
    chr_df$chr_simplified <- rep(NA, length=nrow(chr_df))
    for (gene in chr_df$gene) {
      chr_df[which(chr_df$gene == gene), "chr_name"] <- mapped_genes[which(mapped_genes$Gene == gene), "chromosome_name"]
      chr_df[which(chr_df$gene == gene), "chr_simplified"] <- mapped_genes[which(mapped_genes$Gene == gene), "chr_simplified"]
    }
    chr_df$chr_name <- as.factor(chr_df$chr_name)
    new_chr_list <- append(new_chr_list, list(chr_df))
  }
  names(new_chr_list) <- names(chr_list)
  return(new_chr_list)
}

# 8. Analyze each ct
ProcessCt <- function(main_dir, dis_type, ext, row_col) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  dir.create(out_path, showWarnings = F, recursive = T)
  path <- paste0(main_dir, dis_type, "/outputs/01B_num_DEGs/")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDE(paste(path, sub_ct[ct], sep="/"))
    deg_filt <- list()
    for (k in 1:length(deg)) {
      df_filt <- as.data.frame(rownames(deg[[k]]))
      colnames(df_filt) <- c("Gene")
      deg_filt <- append(deg_filt, list(df_filt))
    }
    names(deg_filt) <- lapply(1:length(names(deg)), function(i) str_replace(names(deg)[i], "_filt", ""))
    for (i in names(deg_filt)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(deg_filt[[i]]))
        names_F <- c(names_F, sub_ct[ct])
      } else {
        df_M <- append(df_M, list(deg_filt[[i]]))
        names_M <- c(names_M, sub_ct[ct])
      }
    }
  }
  names(df_F) <- names_F
  names(df_M) <- names_M
  Annot_F <- lapply(df_F, function(x) Annot.chr.name(x$Gene))
  Annot_M <- lapply(df_M, function(x) Annot.chr.name(x$Gene))
  chr_F <- list()
  chr_M <- list()
  for(i in 1:length(df_F)){
    chr_F_df <- map_chr(df_F[[i]], Annot_F[[i]])
    chr_F <- append(chr_F,list(chr_F_df))
  }
  for(i in 1:length(df_M)){
    chr_M_df <- map_chr(df_M[[i]], Annot_M[[i]])
    chr_M <- append(chr_M,list(chr_M_df))
  }
  names(chr_F) <- names(df_F)
  names(chr_M) <- names(df_M)
  num_chrF <- num_chr_order(chr_F, out_path, "F")
  num_chrM <- num_chr_order(chr_M, out_path, "M")
  #num_chr_list <- list(num_chrF, num_chrM)
  #names(num_chr_list) <- c("F", "M")
  return(list("F" = chr_F, "M" = chr_M))
}

# 9. Plot heatmap of sex-genes across cts
PlotSexHeatmap <- function(main_dir, dis_type, chr_sex, chr, sex, ct_ordered) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  dir.create(out_path, showWarnings = F, recursive = T)
  sexdf <- ExtractSexGenes(chr_sex, chr)
  dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(sexdf$ct))]
  sexdf$ct <- factor(sexdf$ct, dis_ct_ordered)
  sexdf <- sexdf[order(sexdf$ct), ]
  pdf(paste0(out_path, chr, "_genes_heatmap_in_", sex, ".pdf"))
  print(
    ggplot(sexdf, aes(ct, gene)) +
      geom_tile(aes(fill=DEG), color = "light grey") + 
      scale_fill_manual(values = c("y" =  "#F8766D", "n"= "#00BFC4")) +
      scale_color_manual(values = c("y" =  "#F8766D", "n"= "#00BFC4")) +
      labs(x = "Cell types", y = paste0(sex, " genes"), fill = "Expressed", main = sex) +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            axis.text.y = element_text(size=8, colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
    
  )
  dev.off()
}

# 10. Plots heatmaps
PlotSexHmp <- function(main_dir, dis_type, chr_sex_list, ct_ordered) {
  PlotSexHeatmap(main_dir, dis_type, chr_sex_list[[1]], "X", "F", ct_ordered)
  PlotSexHeatmap(main_dir, dis_type, chr_sex_list[[2]], "Y", "M", ct_ordered)
}

# 11. Retrieve p-values from HyperGeom's test done in 02A_HyperGeom
ExtractPval <- function(main_dir, dis_type, sex) {
  in_path <- paste0(main_dir, dis_type, "/outputs/02A_HyperGeom_sex_genes/")
  pval <- read.csv(paste0(in_path, sex, "_HyperGeom_results.csv"))
  names(pval)[names(pval) == 'X.1'] <- 'ct'
  #pval$ct <- as.factor(pval$ct)
  pval$signX <- rep(NA, length(pval$X_enriched_pval))
  pval$signY <- rep(NA, length(pval$Y_enriched_pval))
  pval[which(pval$X_enriched_pval<= 0.05), "signX"] <- "*"
  pval[which(pval$X_enriched_pval<= 0.01), "signX"] <- "**"
  pval[which(pval$X_enriched_pval<= 0.001), "signX"] <- "***"
  pval[which(pval$Y_enriched_pval<= 0.05), "signY"] <- "#"
  pval[which(pval$Y_enriched_pval<= 0.01), "signY"] <- "##"
  pval[which(pval$Y_enriched_pval<= 0.001), "signY"] <- "###"
  return(pval)
}

# 12. Add p-value to df
AddPval <- function(df, pvalX, pvalY) {
  df$pval_HyperGeom <- rep("", length(df$perc))
  for (cells in pvalX$ct) {
    df[which(df$ct==cells & df$chr=="X"),"pval_HyperGeom"] <- pvalX[which(pvalX$ct==cells), "signX"]
  }
  for (cells in pvalY$ct) {
    df[which(df$ct==cells & df$chr=="Y"),"pval_HyperGeom"] <- pvalY[which(pvalY$ct==cells), "signY"]
  }
  return(df)
}

# 13. Plot Heatmap of all DEGs, with hierarchy of chromosome origin
PlotGeneralHeatmap <- function(main_dir, dis_type, chr_sex_list, ct_ordered, group_name) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  dir.create(out_path, showWarnings = F, recursive = T)
  df_sex_list <- ExtractGenes(chr_sex_list) 
  for (sex in names(df_sex_list)) {
    dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(df_sex_list[[sex]]$ct))]
    df_sex_list[[sex]]$ct <- factor(df_sex_list[[sex]]$ct, dis_ct_ordered)
    df_sex_list[[sex]] <- df_sex_list[[sex]][order(df_sex_list[[sex]]$ct), ]
    pdf(paste0(out_path, "all_degs_heatmap_in_", sex, ".pdf"))
    print(
      ggplot(df_sex_list[[sex]], aes(ct, gene)) +
        geom_tile(aes(fill=DEG)) + 
        scale_fill_manual(values = c("y" =  "#F8766D", "n"= "#00BFC4")) +
        scale_color_manual(values = c("y" =  "#F8766D", "n"= "#00BFC4")) +
        labs(x = "Cell types", y = paste0(sex, " DEGs"), fill = "Expressed", title = group_name) +
        facet_wrap(~chr_simplified, scales = "free") +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              plot.title = element_text(size=14, face="bold", colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              axis.text.y = element_blank(),
              axis.ticks.y=element_blank(),
              legend.position = "bottom", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
      
    )
    dev.off()
  }
}

# 14. Plot the fraction of enriched DEGs per chromosome, including or not the HyperGeom p-value
PlotNumChr <- function(main_dir, dis_type, num_chr_genes, ct_ordered, pval_file=F) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  dir.create(out_path, showWarnings = F, recursive = T)
  sexes <- c("F", "M")
  col_palette <- hue_pal()(3)
  for (sex in sexes) {
    if (pval_file) {
      pval_flag <- T
      pval <- ExtractPval(main_dir, dis_type, sex)
      pvalX <- pval[which(pval$X_enriched_pval <= 0.05), ]
      pvalY <- pval[which(pval$Y_enriched_pval <= 0.05), ]
      if (nrow(pvalX) == 0 & nrow(pvalY) == 0) {
        pval_flag <- F
      } 
    } else {
      pvalX <- data.frame()
      pvalX$ct <- vector()
      pvalY <- data.frame()
      pvalY$ct <- vector()
      }
    df <- read.csv(paste0(out_path, sex, "_num_chr.csv"))
    df[,1] <- NULL
    col_factors <- c("ct", "chr")
    df[col_factors] <- lapply(df[col_factors], as.factor) 
    df$perc <- rep(NA, length = length(df$count))
    for (chr in levels(df$chr)) {
      df[which(df$chr==chr), "perc"]  <- df[which(df$chr==chr), "count"] * 100 / num_chr_genes[[chr]]
    }
    if (length(pvalX$ct)>0 | length(pvalY$ct)>0) {
      df <- AddPval(df, pvalX, pvalY)
    }
    dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(df$ct))]
    df$ct <- factor(df$ct, dis_ct_ordered)
    df <- df[order(df$ct), ]
    pdf(paste0(out_path, sex, "_num_chr.pdf"))
    print(
      ggplot(df[which(df$perc>0),], aes(ct, perc, fill=chr)) +
        geom_bar(stat='identity', position=position_dodge2(width = 0.9, preserve = "single")) +
        labs(title=sex, x="", y="Enriched DEGs (%)", fill="Chromosomes") +
        scale_fill_manual(values = c("X" = col_palette[1],
                                     "Y"= col_palette[3],
                                     "Autosome"= col_palette[2])) +
        #{if (length(pvalX$ct)>0) geom_text(size=8, x = pvalX$ct, y = df[which(pvalX$ct %in% df$ct & df$chr=="X"), "perc"] + 0.2, label=pvalX$signX)} +
        #{if (length(pvalY$ct)>0) geom_text(size=8, x = pvalY$ct, y = df[which(pvalY$ct %in% df$ct & df$chr=="Y"), "perc"] + 0.2, label=pvalY$signY)} + 
        {if (pval_flag) geom_text(aes(x = ct, y = perc, label = pval_HyperGeom, fontface = "bold"), nudge_x = 0.15, nudge_y = 0.1)} +
        {if (pval_flag) annotate("text",x= length(levels(df$ct))/2,y=max(df$perc) + 0.2,label="Significance: *: X-chr, #: Y-chr")} +
        {if (pval_flag) coord_cartesian(clip="off") } +
        {if (pval_flag==F) annotate("text", x = length(levels(df$ct))/2, y=max(df$perc), fontface = "bold", label = "NS", colour = "black")} +
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
}

# 15. Saves to CSV output the DEGs from one sex shared by al least 50% of the ct
ExtractSharedGenes <- function(main_dir, dis_type, chr_sex_list, min_shared=0.75) {
  out_path <- paste0(main_dir, dis_type, "/outputs/01C_num_chr/")
  dir.create(out_path, showWarnings = F, recursive = T)
  df_sex_list <- ExtractGenes(chr_sex_list)
  for (sex in names(df_sex_list)) {
    shared_genes <- data.frame()
    thresh <- ceiling(length(levels(df_sex_list[[sex]]$ct)) * min_shared)
    for (gene in unique(df_sex_list[[sex]]$gene)) {
      if (sum(df_sex_list[[sex]][which(df_sex_list[[sex]]$gene==gene), "DEG"]=="y") >= thresh)
        shared_genes <- rbind(shared_genes, df_sex_list[[sex]][which(df_sex_list[[sex]]$gene==gene & df_sex_list[[sex]]$DEG=="y"), c("ct", "gene", "chr_name", "chr_simplified")])
    }
    write.csv(shared_genes, paste0(out_path, "/", sex, "_shared_genes.csv"))
  }
}
