library(tidyr)
library(matrixStats)
library(dplyr)
library(ggplot2)

CalculateSexSD <- function(cell_info_df, top_2000_df) {
  top_list <- list()
  top_names <- vector()
  for (ct_id in unique(cell_info_df$ct)) {
    print(ct_id)
    top_ct_F <- top_2000_df[,c("Genes", cell_info_df[which(cell_info_df$ct==ct_id & cell_info_df$sex=="F"), "cell_id"])]
    top_ct_M <- top_2000_df[,c("Genes", cell_info_df[which(cell_info_df$ct==ct_id & cell_info_df$sex=="M"), "cell_id"])]
    if (ncol(as.data.frame(top_ct_F)) < 3) {
      print(paste0("The cell type ", ct_id, " has less than 2 cells in the F population, therefore there are not enough cells to calculate the SD"))
    } else if (ncol(as.data.frame(top_ct_M)) < 3) {
      print(paste0("The cell type ", ct_id, " has less than 2 in the M population, therefore there are not enough cells to calculate the SD"))
    } else {
      top_ct_F$F_SD <- rowSds(as.matrix(top_ct_F[-1]))
      top_ct_M$M_SD <- rowSds(as.matrix(top_ct_M[-1]))
      top_filt <- merge.data.frame(top_ct_F, top_ct_M, by = "Genes")
      top_filt <- top_filt[,c("Genes", "F_SD", "M_SD")]
      top_list <- append(top_list, list(top_filt))
      top_names <- c(top_names, ct_id)
    }
  }
  names(top_list) <- top_names
  return(top_list)
}

PlotSexSD <- function(main_dir, top_list) {
  plot_path <- paste0(main_dir, "plots/Gene_Variability/")
  dir.create(plot_path, recursive = T, showWarnings = F)
  for (ct in names(top_list)) {
    pdf(paste0(plot_path, ct, "_sex_SD.pdf"))
    print(
      ggplot(top_list[[ct]], aes(F_SD, M_SD)) +
        geom_point() +
        geom_ribbon(aes(F_SD,ymin=0+.9*F_SD, ymax=0+1.1*F_SD), fill="grey70",alpha=.5) +
        geom_abline(intercept = 0, slope = 1, colour="red", linetype = "dashed", linewidth=1) +
        labs(x="Female SD", y="Male SD", title = ct) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              panel.spacing.x=unit(0, "lines"),
              plot.title = element_text(size=12, face="bold", colour = "black"),
              axis.line = element_line(colour = "black"),
              axis.title.x = element_text(size=12, face="bold", colour = "black"),
              axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
              axis.ticks.x=element_blank(),
              axis.title.y = element_text(size=12, face="bold", colour = "black"),
              legend.position = "right", 
              legend.title = element_text(size=12, face="bold", colour = "black"))
    )
    dev.off()
  }
}


SexSD <- function(main_dir, cell_info_df, top_2000_df) {
  top_ct <-CalculateSexSD(cell_info_df, top_2000_df)
  PlotSexSD(main_dir, top_ct)
}