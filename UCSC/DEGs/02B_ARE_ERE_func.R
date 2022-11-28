# 0. Import libraries
library(readxl)
library(stringr)
library(reshape)
library(ggplot2)

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
ReadRawData <- function(main_dir,  sex, ext, row_col) {
  path <- paste0(main_dir,  "/01A_DEGs")
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

ReadRawData2 <- function(main_dir,  sex, ext, row_col) {
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
  return(df_sex)
}

# 3. Calculate ARE df
AREdf <- function(df_sex, sex, ARE_DF) {
  df_ARE <- data.frame(as.factor(names(df_sex)))
  colnames(df_ARE) <- c("ct")
  df_ARE$ct <- sapply(1:length(df_ARE$ct), function(i) str_replace(df_ARE$ct[i], paste0("_", sex, "_filt"), ""))
  df_ARE$bg <- sapply(1:length(names(df_sex)), function(i) nrow(df_sex[[i]]))
  df_ARE$full <- sapply(1:length(names(df_sex)), function(i) length(intersect(rownames(df_sex[[i]]), unique(ARE_DF$fullsites))))
  df_ARE$half <- sapply(1:length(names(df_sex)), function(i) length(intersect(rownames(df_sex[[i]]), unique(ARE_DF$halfsites))))
  df_ARE$hf <- sapply(1:length(names(df_sex)), function(i) length(intersect(rownames(df_sex[[i]]), intersect(unique(ARE_DF$halfsites), unique(ARE_DF$fullsites)))))
  df_ARE$full <- df_ARE$full - df_ARE$hf
  df_ARE$half <- df_ARE$half - df_ARE$hf
  df_ARE$no_overlap <- df_ARE$bg - df_ARE$full - df_ARE$half - df_ARE$hf
  return(df_ARE)
}

# 4. Calculate percentages of ARE sites
AREdfPerc <- function(main_dir,  df_ARE, sex) {
  df_ARE <- transform(df_ARE, full_perc = full * 100 / bg)
  df_ARE <- transform(df_ARE, half_perc = half * 100 / bg)
  df_ARE <- transform(df_ARE, hf_perc = hf * 100 / bg)
  df_ARE <- transform(df_ARE, no_overlap_perc = no_overlap * 100 / bg)
  dir.create(paste(main_dir,  "02B_ARE_ERE", sep="/"), showWarnings = FALSE)
  write.csv(df_ARE, paste0(main_dir, "/02B_ARE_ERE/", sex, "_ARE_sites.csv"))
  df_ARE_perc <- df_ARE[, c(1, 7:10)]
  df_ARE_perc <- melt(df_ARE_perc, id.vars = "ct")
  names(df_ARE_perc)[names(df_ARE_perc) == 'value'] <- 'percent'
  names(df_ARE_perc)[names(df_ARE_perc) == 'variable'] <- 'sites'
  col_factors <- c("ct")
  df_ARE_perc[col_factors] <- lapply(df_ARE_perc[col_factors], as.factor) 
  levels(df_ARE_perc$sites) <- c('Full', 'Half', 'Half-Full', 'None')
  return(df_ARE_perc)
}

# 5. Calculate ERE df
EREdf <- function(df_sex, sex, ERE_gene) {
  df_ERE <- data.frame(as.factor(names(df_sex)))
  colnames(df_ERE) <- c("ct")
  df_ERE$ct <- sapply(1:length(df_ERE$ct), function(i) str_replace(df_ERE$ct[i], paste0("_", sex, "_filt"), ""))
  df_ERE$bg <- sapply(1:length(names(df_sex)), function(i) nrow(df_sex[[i]]))
  df_ERE$ERE_overlap <- sapply(1:length(names(df_sex)), function(i) length(intersect(rownames(df_sex[[i]]), unique(ERE_gene))))
  df_ERE$no_overlap <- df_ERE$bg - df_ERE$ERE_overlap
  return(df_ERE)
}

# 6. Calculate percentages of ARE sites
EREdfPerc <- function(main_dir,  df_ERE, sex) {
  df_ERE <- transform(df_ERE, ERE_perc = ERE_overlap * 100 / bg)
  df_ERE <- transform(df_ERE, no_overlap_perc = no_overlap * 100 / bg)
  dir.create(paste(main_dir,  "02B_ARE_ERE", sep="/"), showWarnings = FALSE)
  write.csv(df_ERE, paste0(main_dir,   "/02B_ARE_ERE/", sex, "_ERE_sites.csv"))
  df_ERE_perc <- df_ERE[, c(1, 5:6)]
  df_ERE_perc <- melt(df_ERE_perc, id.vars = "ct")
  names(df_ERE_perc)[names(df_ERE_perc) == 'value'] <- 'percent'
  names(df_ERE_perc)[names(df_ERE_perc) == 'variable'] <- 'sites'
  col_factors <- c("ct")
  df_ERE_perc[col_factors] <- lapply(df_ERE_perc[col_factors], as.factor) 
  levels(df_ERE_perc$sites) <- c("ERE", "None")
  return(df_ERE_perc)
}

# 5. Plot ARE or ERE sites
PlotAREERE <- function(main_dir,  df_sex, sex, are_ere, ct_ordered) {
  if (are_ere=="ARE") {
    col_palette <- c("#39B600", "#9590FF","#D376FF" , "#FD61D1")
  } else if (are_ere=="ERE") {
    col_palette <- c("#39B600", "#9590FF")
  }
  dis_ct_ordered <- ct_ordered[which(ct_ordered %in% levels(df_sex$ct))]
  df_sex$ct <- factor(df_sex$ct, dis_ct_ordered)
  df_sex <- df_sex[order(df_sex$ct), ]
  pdf(paste0(main_dir,   "/02B_ARE_ERE/", sex, "_", are_ere, "_sites.pdf"))
  print(ggplot(df_sex, aes(ct, percent, fill=sites)) +
          geom_bar(stat="identity", position="stack", color="black") + 
          labs(title=sex, x="Cell types", y=paste0("% of ", are_ere, " sites"), fill=paste0("Overlap ", are_ere , " sites")) +
          {if (are_ere=="ARE") scale_fill_manual(values = c('Full' = col_palette[4] , 'Half' = col_palette[3], 'Half-Full' = col_palette[2], 'None' = col_palette[1]))} +
          {if (are_ere=="ERE") scale_fill_manual(values = c("ERE" = col_palette[2], "None" = col_palette[1]))} +
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


# 7. MAIN
AnalysisARE_ERE <- function(main_dir,  ARE_DF, ERE_gene, ct_order) {
  sexes <- c("F", "M")
  for (sex in sexes) {
    df_sex <- ReadRawData2(main_dir,  sex)
    ARE_sex <- AREdf(df_sex, sex, ARE_DF)
    ARE_sex <- AREdfPerc(main_dir,  ARE_sex, sex)
    ERE_sex <- EREdf(df_sex, sex, ERE_gene)
    ERE_sex <- EREdfPerc(main_dir,  ERE_sex, sex)
    PlotAREERE(main_dir,  ARE_sex, sex, "ARE", ct_order)
    PlotAREERE(main_dir,  ERE_sex, sex, "ERE", ct_order)
  }
}