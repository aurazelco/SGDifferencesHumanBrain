library(ggplot2)
library(readxl)
library(stringr)
library(tidyr)

SCENICct <- function(main_dir, dis_type) {
  ct_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells/")
  ct_files <- list.files(path = ct_path, pattern = "\\.tsv$",full.names = T)
  cts <- lapply(ct_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(cts) <- list.files(path = ct_path, pattern = "\\.tsv$",full.names = F)
  names(cts) <- sapply(1:length(names(cts)), function(x) str_remove(names(cts)[x], ".tsv"))
  for (df in names(cts)) {
    cts[[df]]$TF_TG <- paste(cts[[df]]$TF, cts[[df]]$target, sep="_")
    cts[[df]]$TF_TG <- as.factor(cts[[df]]$TF_TG)
  }
  cts_sort <- lapply(1:length(names(cts)), function(x) cts[[x]][order(-cts[[x]]$importance),] )
  names(cts_sort) <- names(cts)
  return(cts_sort)
}

SCENICall <- function(main_dir, dis_type) {
  all_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells_all/")
  all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
  alls <- lapply(all_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(alls) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
  names(alls) <- sapply(1:length(names(alls)), function(x) str_remove(names(alls)[x], ".tsv"))
  for (df in names(alls)) {
    alls[[df]]$TF_TG <- paste(alls[[df]]$TF, alls[[df]]$target, sep="_")
    alls[[df]]$TF_TG <- as.factor(alls[[df]]$TF_TG)
  }
  alls_sort <- lapply(1:length(names(alls)), function(x) alls[[x]][order(-alls[[x]]$importance),] )
  names(alls_sort) <-names(alls) 
  return(alls_sort)
}

SCENICoverlap <- function(main_dir, dis_type) {
  cts_sort <- SCENICct(main_dir, dis_type)
  alls_sort <- SCENICall(main_dir, dis_type)
  return(list(cts_sort, alls_sort))
}

PlotOverlapRuns <- function(main_dir, dis_type, sort_list, top_list, ct_list) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  cts_sort <- sort_list[[1]]
  alls_sort <- sort_list[[2]]
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  for (k in top_list) {
    print(k)
    if (k!="no") {
      k <- as.numeric(k)
      ct <- vector()
      tot <- vector()
      common <- vector()
      for (id in ct_list) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        ct <- c(ct, id)
        common <- c(common, length(intersect(sub_list[[1]][1:k, "TF_TG"], intersect(sub_list[[2]][1:k, "TF_TG"], sub_list[[3]][1:k, "TF_TG"]))))
      }
      tot <- rep(k, length(ct))
      overlap_ct <- data.frame(ct, tot, common)
      overlap_ct$perc <- overlap_ct$common*100 / overlap_ct$tot
      overlap_ct$ct <- as.factor(overlap_ct$ct)
      common_all <- length(intersect(alls_sort[[1]][1:k, "TF_TG"], intersect(alls_sort[[2]][1:k, "TF_TG"], alls_sort[[3]][1:k, "TF_TG"])))
      overlap_alls <- data.frame(c("All"), k, common_all, common_all*100 / as.numeric(k))
      colnames(overlap_alls) <- colnames(overlap_ct)
      overlap <- rbind(overlap_alls, overlap_ct)
      overlap$ct <- as.factor(overlap$ct)
      pdf(paste0(main_dir, dis_type, "/3_plots/overlap_from_GRN_all_", k, "_TF_TG_pairs.pdf"))
      print(
        ggplot(overlap, aes(ct, perc, fill=ct)) +
          geom_bar(stat="identity") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    } else {
      ct <- vector()
      tot <- vector()
      common <- vector()
      for (id in ct_list) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        ct <- c(ct, id)
        tot <- c(tot, length(unique(sub_list[[1]][, "TF_TG"], sub_list[[2]][, "TF_TG"], sub_list[[3]][, "TF_TG"])))
        common <- c(common, length(intersect(sub_list[[1]][, "TF_TG"], intersect(sub_list[[2]][, "TF_TG"], sub_list[[3]][, "TF_TG"]))))
      }
      overlap_ct <- data.frame(ct, tot, common)
      overlap_ct$perc <- overlap_ct$common*100 / overlap_ct$tot
      overlap_ct$ct <- as.factor(overlap_ct$ct)
      tot_all <- length(unique(alls_sort[[1]][, "TF_TG"], alls_sort[[2]][, "TF_TG"], alls_sort[[3]][, "TF_TG"]))
      common_all <- length(intersect(alls_sort[[1]][, "TF_TG"], intersect(alls_sort[[2]][, "TF_TG"], alls_sort[[3]][, "TF_TG"])))
      overlap_alls <- data.frame(c("All"), tot_all, common_all, common_all*100 / tot_all)
      colnames(overlap_alls) <- colnames(overlap_ct)
      overlap <- rbind(overlap_alls, overlap_ct)
      overlap$ct <- as.factor(overlap$ct)
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_all_TF_TG_pairs.pdf"))
      print(
        ggplot(overlap, aes(ct, perc, fill=ct)) +
          geom_bar(stat="identity") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
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
}

PlotscGRNomOverlap <- function(main_dir, dis_type, sort_list, top_scGRNom, ct_scGRNom, flag) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  cts_sort <- sort_list[[1]]
  alls_sort <- sort_list[[2]]
  for (f in flag) {
    if (f == "no") {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (id in names(ct_scGRNom)) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        names(sub_list) <- sapply(1:length(names(sub_list)), function(x) str_replace(names(sub_list)[x], id, ct_scGRNom[[id]]))
        sub_scGRNom_names <- names(scGRNom)[which(grepl(ct_scGRNom[[id]], names(scGRNom)))]
        sub_scGRNom <- scGRNom[sub_scGRNom_names]
        for (scGRNom_id in names(sub_scGRNom)) {
          for (id_sub in names(sub_list)) {
            sort_sub <- sub_list[[id_sub]][order(-sub_list[[id_sub]]$importance),]
            sort_sub <- sort_sub[1:top_scGRNom[[scGRNom_id]], ]
            scGRNom_groups <- c(scGRNom_groups, paste(id_sub, scGRNom_id, sep="_"))
            tot_scGRNom <- c(tot_scGRNom, length(unique(sort_sub$TF_TG, sub_scGRNom[scGRNom_id]$TF_TG)))
            common_scGRNom <- c(common_scGRNom, length(intersect(sort_sub$TF_TG, sub_scGRNom[[scGRNom_id]]$TF_TG)))
          }
        }
      }
      overlap_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_scGRNom$perc <- overlap_scGRNom$common_scGRNom * 100 / overlap_scGRNom$tot_scGRNom
      overlap_scGRNom$scGRNom_groups <- str_replace_all(overlap_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_scGRNom <- separate(overlap_scGRNom, scGRNom_groups, into=c("disco", "run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("disco", "run", "GRN_ct", "openchrom")
      overlap_scGRNom[col_factors] <- lapply(overlap_scGRNom[col_factors], as.factor)
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black")))
      dev.off()
    } else {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (id in names(ct_scGRNom)) {
        sub_names <- names(cts_sort)[which(grepl(id, names(cts_sort)))]
        sub_list <- cts_sort[sub_names]
        names(sub_list) <- sapply(1:length(names(sub_list)), function(x) str_replace(names(sub_list)[x], id, ct_scGRNom[[id]]))
        sub_scGRNom_names <- names(scGRNom)[which(grepl(ct_scGRNom[[id]], names(scGRNom)))]
        sub_scGRNom <- scGRNom[sub_scGRNom_names]
        for (scGRNom_id in names(sub_scGRNom)) {
          for (id_sub in names(sub_list)) {
            sort_sub <- sub_list[[id_sub]][order(-sub_list[[id_sub]]$importance),]
            sort_sub <- sort_sub[1:top_scGRNom[[scGRNom_id]], ]
            scGRNom_groups <- c(scGRNom_groups, paste(id_sub, scGRNom_id, sep="_"))
            tot_scGRNom <- c(tot_scGRNom, top_scGRNom[[scGRNom_id]])
            common_scGRNom <- c(common_scGRNom, length(intersect(sort_sub$TF_TG, sub_scGRNom[[scGRNom_id]]$TF_TG)))
          }
        }
      }
      overlap_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_scGRNom$perc <- overlap_scGRNom$common_scGRNom * 100 / overlap_scGRNom$tot_scGRNom
      overlap_scGRNom$scGRNom_groups <- str_replace_all(overlap_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_scGRNom <- separate(overlap_scGRNom, scGRNom_groups, into=c("disco", "run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("disco", "run", "GRN_ct", "openchrom")
      overlap_scGRNom[col_factors] <- lapply(overlap_scGRNom[col_factors], as.factor)
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overlap_from_GRN_cts_scGRNom.pdf"))
      print(
        ggplot(overlap_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black")))
      dev.off()
    }
  }
  for (f in flag) {
    if (f=="no") {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (a in names(alls_sort)) {
        for (c in names(scGRNom)) {
          scGRNom_groups <- c(scGRNom_groups, paste(a, c, sep="_"))
          a_sort <- alls_sort[[a]][1:top_scGRNom[[c]], "TF_TG"]
          common_scGRNom <- c(common_scGRNom, length(intersect(a_sort, sub_scGRNom[[c]]$TF_TG)))
          tot_scGRNom <- c(tot_scGRNom, length(unique(a_sort, sub_scGRNom[[c]]$TF_TG)))
        }
      }
      overlap_alls_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_alls_scGRNom$scGRNom_groups <- str_replace_all(overlap_alls_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_alls_scGRNom <- separate(overlap_alls_scGRNom, scGRNom_groups, into=c("run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("run", "GRN_ct", "openchrom")
      overlap_alls_scGRNom[col_factors] <- lapply(overlap_alls_scGRNom[col_factors], as.factor)
      overlap_alls_scGRNom$perc <- overlap_alls_scGRNom$common_scGRNom * 100 / overlap_alls_scGRNom$tot_scGRNom
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overall_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overall among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
    } else {
      scGRNom_groups <- vector()
      tot_scGRNom <- vector()
      common_scGRNom <- vector()
      for (a in names(alls_sort)) {
        for (c in names(scGRNom)) {
          scGRNom_groups <- c(scGRNom_groups, paste(a, c, sep="_"))
          tot_scGRNom <- c(tot_scGRNom, top_scGRNom[[c]])
          a_sort <- alls_sort[[a]][1:top_scGRNom[[c]], "TF_TG"]
          common_scGRNom <- c(common_scGRNom, length(intersect(a_sort, sub_scGRNom[[c]]$TF_TG)))
        }
      }
      overlap_alls_scGRNom <- data.frame(scGRNom_groups, tot_scGRNom, common_scGRNom)
      overlap_alls_scGRNom$scGRNom_groups <- str_replace_all(overlap_alls_scGRNom$scGRNom_groups, c("_Mic"="/Mic", "_Oli"="/Oli", "_with"="/with"))
      overlap_alls_scGRNom <- separate(overlap_alls_scGRNom, scGRNom_groups, into=c("run", "GRN_ct", "openchrom"), remove=F, sep="/")
      col_factors <- c("run", "GRN_ct", "openchrom")
      overlap_alls_scGRNom[col_factors] <- lapply(overlap_alls_scGRNom[col_factors], as.factor)
      overlap_alls_scGRNom$perc <- overlap_alls_scGRNom$common_scGRNom * 100 / overlap_alls_scGRNom$tot_scGRNom
      pdf(paste0(main_dir, dis_type, "/3_plots/unsorted_overlap_from_GRN_alls_scGRNom.pdf"))
      print(
        ggplot(overlap_alls_scGRNom, aes(run, perc, fill=openchrom)) +
          geom_bar(stat="identity", position="dodge", color="black") +
          labs(x="Cell types with num cell > 500", y="% Overlap among runs and scGRNom analysis", fill="Groups", title="GSE157827_F_Normal") +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
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
}

SCENICresults <- function(main_dir, dis_type) {
  dir.create(paste0(main_dir, dis_type, "/3_plots"), showWarnings = F)
  all_path <- paste0(main_dir, dis_type, "/1_GRN/100_sampled_cells_all/")
  all_files <- list.files(path = all_path, pattern = "\\.tsv$",full.names = T)
  all <- lapply(all_files,function(x) {
    read.table(file = x, 
               sep = '\t', 
               header = TRUE)
  })
  names(all) <- list.files(path = all_path, pattern = "\\.tsv$",full.names = F)
  only_1 <- names(all)[which(grepl("_1", names(all)))]
  all_1 <- all[only_1]
  names(all_1) <- str_remove_all(names(all_1), "_1.tsv")
  all_sort <- lapply(1:length(names(all_1)), function(x) all_1[[x]][order(-all_1[[x]]$importance),])
  names(all_sort) <- names(all_1)
  return(all_sort)
}

SCENICInput <- function(main_dir, dis_type, cell_info_df) {
  input_dfs_path <-  paste0(main_dir, dis_type, "/0_input_dfs/sampled_100_cells_all")
  input_dfs_files <- list.files(path = input_dfs_path, full.names = T)
  input_dfs <- lapply(input_dfs_files,function(x) {
    read.table(file = x, 
               sep = ',', 
               header = TRUE)
  })
  names(input_dfs) <- list.files(path = input_dfs_path, full.names = F)
  names(input_dfs) <- str_remove_all(names(input_dfs), ".csv")
  only_1 <- names(input_dfs)[which(grepl("_1", names(input_dfs)))]
  input_dfs <- input_dfs[only_1]
  names(input_dfs) <- str_remove_all(names(input_dfs), "_1")
  for (id in names(input_dfs)) {
    colnames(input_dfs[[id]]) <- str_replace_all(colnames(input_dfs[[id]]), "[.]", "-")
    sub_info <- subset(cell_info_df, cell_id %in% colnames(input_dfs[[id]])[-1])
    ct_info <- sub_info$ct
    colnames(input_dfs[[id]]) <- c("Genes", ct_info)
  }
  return(input_dfs)
}

PlotTfTg <- function(main_dir, dis_type, all_sort, input_dfs, ct_ordered, top_TF){
  for (id in names(all_sort)) {
    for (k in top_TF) {
      tf_id <- all_sort[[id]][1:k, "TF"]
      sub_tf <- subset(input_dfs[[id]], Genes %in% tf_id)
      sub_tf <-  reshape2::melt(sub_tf)
      colnames(sub_tf) <- c("TF", "ct", "value")
      sub_tf$ct <- as.character(sub_tf$ct)
      sub_tf$ct <- str_remove_all(sub_tf$ct, "[.][:digit:]+")
      sub_tf[c("TF", "ct")] <- lapply(sub_tf[c("TF", "ct")], as.factor)
      sub_tf_order <- ct_ordered[which(ct_ordered %in% levels(sub_tf$ct))]
      sub_tf$ct <- factor(sub_tf$ct, sub_tf_order)
      sub_tf <- sub_tf[order(sub_tf$ct), ]
      pdf(paste0(main_dir, dis_type, "/3_plots/TFs_hmp_",k, "_top_", id, ".pdf"))
      print(
        ggplot(sub_tf, aes(ct, TF, fill=value)) +
          geom_tile() + 
          coord_fixed() +
          scale_fill_gradient(low="blue", high="red") +
          labs(y="TFs", fill="RNA Expression", title=id) +
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                plot.title = element_text(size=12, face="bold", colour = "black"),
                axis.line = element_line(colour = "black"),
                axis.title.x = element_blank(),
                axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
                axis.ticks.x=element_blank(),
                axis.title.y = element_text(size=12, face="bold", colour = "black"),
                legend.position = "bottom", 
                legend.title = element_text(size=12, face="bold", colour = "black"))
      )
      dev.off()
      tg_id <- all_sort[[id]][1:k, "target"]
      sub_tg <- subset(input_dfs[[id]], Genes %in% tg_id)
      if (nrow(sub_tg) > 0) {
        sub_tg <-  reshape2::melt(sub_tg)
        colnames(sub_tg) <- c("TG", "ct", "value")
        sub_tg$ct <- as.character(sub_tg$ct)
        sub_tg$ct <- str_remove_all(sub_tg$ct, "[.][:digit:]+")
        sub_tg[c("TG", "ct")] <- lapply(sub_tg[c("TG", "ct")], as.factor)
        sub_tg_order <- ct_ordered[which(ct_ordered %in% levels(sub_tg$ct))]
        sub_tg$ct <- factor(sub_tg$ct, sub_tg_order)
        sub_tg <- sub_tg[order(sub_tg$ct), ]
        if (k==100) {
          pdf(paste0(main_dir, dis_type, "/3_plots/TGs_hmp_",k, "_top_", id, ".pdf"),
              height=11)
        } else {
          pdf(paste0(main_dir, dis_type, "/3_plots/TGs_hmp_",k, "_top_", id, ".pdf"),
              height=15)
        }
        print(
          ggplot(sub_tg, aes(ct, TG, fill=value)) +
            geom_tile() + 
            coord_fixed() +
            scale_fill_gradient(low="blue", high="red") +
            labs(y="TGs", fill="RNA Expression", title=id) +
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(), 
                  plot.title = element_text(size=12, face="bold", colour = "black"),
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
  }
}
