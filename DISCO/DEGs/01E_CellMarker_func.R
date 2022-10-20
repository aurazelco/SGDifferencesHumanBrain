################ USER-DEFINED FUNCTIONS

# 0. Import libraries
library(ggplot2)
#library(dplyr)
library(stringr)

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

# 2. Extract cell markers for each ct class
CellMarkers <- function(main_dir, cm_ct_list) {
  cm1 <- read.csv(paste0(main_dir, "CellMarker.csv"))
  cm2 <- read.csv(paste0(main_dir, "CellMarker-2.csv"))
  cm3 <- read.csv(paste0(main_dir, "CellMarker-3.csv"))
  cm_df <- rbind(cm1, cm2, cm3)
  rm(cm1, cm2, cm3)
  cm_df <- cm_df[, -c(6:8)]
  col_factors <- c("Species",
                   "Tissue",
                   "Cell.Type",
                   "Cancer")
  cm_df[col_factors] <- lapply(cm_df[col_factors], as.factor) 
  cm_df <- subset(cm_df, subset = (Cancer == "Normal"))
  cm_df <- subset(cm_df, subset = (Tissue == c("Brain") | Tissue == c("Dorsolateral prefrontal cortex")))
  cm_df$Cancer <- droplevels(cm_df$Cancer)
  cm_df$Tissue <- droplevels(cm_df$Tissue)
  cm_df$Cell.Type <- droplevels(cm_df$Cell.Type)
  cm_df$ct <- rep(NA, nrow(cm_df))
  
  for (ct in levels(cm_df$Cell.Type)) {
    cm_df[which(cm_df$Cell.Type==ct), "ct"] <- cm_ct_list[[ct]]
  }
  cm_df$ct <- as.factor(cm_df$ct)
  markers <- list()
  for (i in levels(cm_df$ct)) {
    markers <- append(markers, list((cm_df[which(cm_df$ct==i),"Cell.Marker"])))
  }
  names(markers) <- levels(cm_df$ct)
  for (ct in names(markers)) {
    gene_list <- vector()
    for (i in 1:length(markers[[ct]])) {
      gene_list <- c(gene_list, str_split(markers[[ct]][i], ", "))
    }
    markers[[ct]] <- unique(unlist(gene_list))
  }
  return(markers)
}

# 3. Calculates the percentage of DEGs which are also markers for the cts
CalcPercMarkers <- function(df_sex, markers, markers_ct_names, ct_order) {
  if (length(names(df_sex[[1]])) == length(names(df_sex[[2]]))) {
    df_markers <- as.data.frame(rep(unlist(names(df_sex[[1]])), 2))
    colnames(df_markers) <- c("ct")
    df_markers$sex <- rep(names(df_sex), each = length(names(df_sex[[1]])))
  } else {
    ct_names <- unlist(names(df_sex[[1]]), names(df_sex[[2]])) 
    df_markers <- as.data.frame(ct_names)
    colnames(df_markers) <- c("ct")
    df_markers$sex <- c(rep(names(df_sex)[1], length(names(df_sex[[1]]))), rep(names(df_sex)[2], length(names(df_sex[[2]]))))
  }
  df_markers$markers_count <- rep(NA, nrow(df_markers))
  df_markers$markers_perc <- rep(NA, nrow(df_markers))
  for (sex in names(df_sex)) {
    for (ct_name in names(df_sex[[sex]])) {
      df_markers[which(df_markers$ct==ct_name & df_markers$sex==sex),"markers_count"] <- length(intersect(df_sex[[sex]][[ct_name]][[sex]], markers[[data_ct[[ct_name]]]]))
      df_markers[which(df_markers$ct==ct_name & df_markers$sex==sex),"markers_perc"] <- length(intersect(df_sex[[sex]][[ct_name]][[sex]], markers[[data_ct[[ct_name]]]]))*100/length(markers[[data_ct[[ct_name]]]])
    }
  }
  col_factors <- c("ct", "sex")
  df_markers[col_factors] <- lapply(df_markers[col_factors], as.factor)
  dis_ct_order <- ct_order[which(ct_order %in% levels(df_markers$ct))]
  df_markers$ct <- factor(df_markers$ct, dis_ct_order)
  df_markers <- df_markers[order(df_markers$ct), ]
  return(df_markers)
}

# 4. Plots the results
PlotMarkers <- function(main_dir, dis_type, df_markers) {
  dir.create(paste0(main_dir, dis_type, "/01E_CellMarker"), showWarnings = FALSE)
  pdf(paste0(main_dir, dis_type, "/01E_CellMarker/CellMarker_count.pdf"))
  print(
    ggplot(df_markers[which(df_markers$markers_count>0), ], aes(ct, markers_count, fill=sex)) +
      geom_bar(stat='identity', position=position_dodge2(width = 0.9, preserve = "single")) +
      labs(x="Cell types", y="Abs counts of Markers", fill="Sex") +
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
  pdf(paste0(main_dir, dis_type, "/01E_CellMarker/CellMarker_percentage.pdf"))
  print(
  ggplot(df_markers[which(df_markers$markers_perc>0), ], aes(ct, markers_perc, fill=sex)) +
    geom_bar(stat='identity', position=position_dodge2(width = 0.9, preserve = "single")) +
    labs(x="Cell types", y="% Markers", fill="Sex") +
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

# 5. Combines prevoius functions into one
PlotCMresults <- function(main_dir, dis_type, cm_dir, cm_ct_list, markers_ct_names, ct_order, row_col) {
  path <- paste0(main_dir, dis_type, "/01B_num_DEGs")
  sub_ct <- list.dirs(path, recursive=FALSE, full.names = FALSE)
  df_F <- list()
  df_M <- list()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportSignDE(paste(path, sub_ct[ct], sep="/"))
    for (i in names(deg)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(deg[[i]]))
      } else {
        df_M <- append(df_M, list(deg[[i]]))
      }
    }
  }
  names(df_F) <- sub_ct
  names(df_M) <- sub_ct
  df_sex <- (list("F" = df_F, "M" = df_M))
  markers <- CellMarkers(cm_dir, cm_ct_list)
  df_markers <- CalcPercMarkers(df_sex, markers, markers_ct_names, ct_order)
  PlotMarkers(main_dir, dis_type, df_markers)
}