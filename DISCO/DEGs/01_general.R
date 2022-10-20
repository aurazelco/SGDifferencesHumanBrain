main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

library(ggplot2)

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

num_filt <- read.csv(paste0(main, "final_filt.csv"))
num_filt[,1] <- NULL

col_factors <- c("proj", "sex", "disease", "ct", "og", "idents")
num_filt[col_factors] <- lapply(num_filt[col_factors], as.factor)  

num_filt$ct <- factor(num_filt$ct, levels=c(
  "L2/3 EN",               
  "L4 EN",   
  "PLCH1 L4/5 EN", 
  "TSHZ2 L4/5 EN", 
  "L5 EN",       
  "L5/6 EN",       
  "L5b EN",     
  "L6 EN",     
  "pyramidal neuron", 
  "CXCL14 IN",  
  "PVALB IN",                    
  "SST IN",
  "SV2C IN",               
  "VIP IN",  
  "EC", 
  "fibrous astrocyte",
  "protoplasmic astrocyte",
  "OPC", 
  "oligodendrocyte",           
  "microglia",   
  "T"
  ))


# Plot number of cells per sex and ct
for (dis_type in sub_disease) {
  filt_df <- subset(num_filt, disease == dis_type)
  pdf(paste0(main, dis_type, "/num_ct_sex.pdf"))
  print(
    ggplot(filt_df, aes(ct, count, fill=ct)) +
      geom_bar(stat="identity") +
      facet_wrap(~sex, scales="free") +
      labs(y="Number of cells", fill="") +
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
  pdf(paste0(main, dis_type, "/num_ct_sex_proj.pdf"))
  print(
    ggplot(filt_df, aes(ct, count, fill=ct)) +
      geom_bar(stat="identity") +
      facet_wrap(~sex*proj, scales="free") +
      labs(y="Number of cells", fill="") +
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
