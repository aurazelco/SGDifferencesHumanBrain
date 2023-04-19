library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(ggpubr)


out_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/suppl_files/"
plot_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Extra_figures/"


disco <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/num_proj_sex_ct.csv")

velm <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_sex_ct_per_age.csv")
#eze_nowa <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Eze_Nowakowski_num_cells.csv")

disco$X <- NULL
disco$disease <- str_replace_all(disco$disease, "Normal", "Healthy")
velm$X <- NULL
#eze_nowa$X <- NULL
#eze_nowa <- cbind(rep("Eze_Nowakowski_2nd_trimester", nrow(eze_nowa)), eze_nowa)

#colnames(eze_nowa) <- c("proj", "sex", "ct", "count")
colnames(velm) <- c("proj", "sex", "ct", "count")
#eze_nowa$disease <- rep("Normal", nrow(eze_nowa))
velm$disease <- rep("Healthy", nrow(velm))

#eze_nowa <- eze_nowa %>% relocate(disease, .after = sex)
velm <- velm %>% relocate(disease, .after = sex)

velm$proj <- paste("Velmeshev_2022", velm$proj, sep = "_")
disco$proj <- paste("DISCO", disco$proj, sep = "_")


#all_ds_ct <- rbind(eze_nowa, velm, disco)
all_ds_ct <- rbind(velm, disco)
all_ds_ct$sex <- str_replace_all(all_ds_ct$sex, c("Female"="F", "Male"="M"))

sex_count <- vector()
id <- vector()
for (ds in unique(all_ds_ct$proj)) {
  for (dis in unique(all_ds_ct[which(all_ds_ct$proj==ds), "disease"])) {
    sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds & all_ds_ct$disease==dis & all_ds_ct$sex=="F"), "count"]))
    sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds & all_ds_ct$disease==dis &  all_ds_ct$sex=="M"), "count"]))
    id <- c(id, paste(ds, dis, "F", sep = "/"), paste(ds, dis, "M", sep = "/"))
  }
}

all_ds <- data.frame("id"=id, "count"=sex_count)
all_ds <- separate(all_ds, id, into = c("proj","disease", "sex"), sep = "/")

write.csv(all_ds, paste0(out_path, "num_cells.csv"))

order_proj <- c("Eze_Nowakowski_2nd_trimester", "Velmeshev_2022_2nd trimester",
                "Velmeshev_2022_3rd trimester", "Velmeshev_2022_0-1 years", "Velmeshev_2022_1-2 years",
                "Velmeshev_2022_2-4 years", "Velmeshev_2022_10-20 years", "Velmeshev_2022_Adult",
                "DISCO_GSE157827", "DISCO_GSE174367","DISCO_PRJNA544731")

norm <- ggplot(all_ds[which(all_ds$disease=="Healthy"),], aes(factor(proj, order_proj[-1]), count, fill=sex)) +
  geom_bar(stat = "identity", color="black", position = "dodge") +
  labs(x="Datasets", y="Number of cells", fill="Sex") +
  facet_wrap(~factor(disease, levels=c("Healthy", "Alzheimer's disease", "Multiple Sclerosis")), 
             scales = "free_x") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=12, face="bold", colour = "black"),
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"))

dis <- ggplot(all_ds[which(all_ds$disease!="Healthy"),], aes(factor(proj, order_proj), count, fill=sex)) +
  geom_bar(stat = "identity", color="black", position = "dodge") +
  labs(x="Datasets", y="Number of cells", fill="Sex") +
  facet_wrap(~factor(disease, levels=c("Healthy", "Alzheimer's disease", "Multiple Sclerosis")), 
             scales = "free_x") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        strip.text = element_text(size=12, face="bold", colour = "black"),
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=8, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"))


fig <- ggarrange(norm, dis, nrow = 2, common.legend = T, legend = "bottom")

png(paste0(plot_path, "num_cells.png"), height = 15, width = 10, units="in", res = 1200 )
print(fig)
dev.off()

unified_annotation <- c("CXCL14 IN" = "Interneurons",
                        "EC" = "Endothelial cells",
                        "fibrous astrocyte"  = "Astrocytes",
                        "L2_3 EN" = "Excitatory neurons", 
                        "L4 EN" = "Excitatory neurons",
                        "L5 EN" = "Excitatory neurons",
                        "L5_6 EN" = "Excitatory neurons",
                        "L5b EN" = "Excitatory neurons",
                        "L6 EN" = "Excitatory neurons",                
                        "microglia" = "Microglia", 
                        "Oligodendrocyte" =  "Oligodendrocytes",      
                        "OPC" = "OPCs",                  
                        "PLCH1 L4_5 EN" = "Excitatory neurons", 
                        "protoplasmic astrocyte" = "Astrocytes",
                        "PVALB IN"  = "Interneurons",            
                        "pyramidal neuron"  = "Excitatory neurons",
                        "SST IN" = "Interneurons",   
                        "SV2C IN"  = "Interneurons",   
                        "T" = "T cells",
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
                        "IPC"="IPCs",
                        "Mesenchymal" = "Mesenchymal",      
                        "Neuroepithelial" =     "Neuroepithelial",
                        "Neuronal" = "Neurons",            
                        "Other"    = "Other",                
                        "Radial Glial"     = "Radial Glia",       
                        "Astrocytes" = "Astrocytes",        
                        "Excitatory neurons"  = "Excitatory neurons",
                        "Interneurons"   = "Interneurons",     
                        "Microglia"  = "Microglia",         
                        "Oligodendrocytes" = "Oligodendrocytes",
                        "OPCs" = "OPCs",            
                        "Unknown" = "Unknown",           
                        "Vascular cells" = "Vascular cells",     
                        "Dorsal progenitors"  = "Dorsal progenitors" ,   
                        "Ventral progenitors" = "Ventral progenitors")
names(unified_annotation) <- tolower(names(unified_annotation))

all_ds_ct$ct <- str_replace_all(all_ds_ct$ct, "/", "_")

all_ds_ct$ct_common <- rep(NA, nrow(all_ds_ct))
for (ct in tolower(unique(all_ds_ct$ct))) {
  all_ds_ct[which(tolower(all_ds_ct$ct)==ct), "ct_common"] <- unified_annotation[ct]
}

sex_count <- vector()
id <- vector()
for (ds in unique(all_ds_ct$proj)) {
  for (dis in unique(all_ds_ct[which(all_ds_ct$proj==ds), "disease"])) {
    for (ct in unique(all_ds_ct[which(all_ds_ct$proj==ds & all_ds_ct$disease==dis), "ct_common"])) {
      sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds & all_ds_ct$disease==dis & all_ds_ct$ct_common==ct & all_ds_ct$sex=="F"), "count"]))
      sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds & all_ds_ct$disease==dis & all_ds_ct$ct_common==ct &  all_ds_ct$sex=="M"), "count"]))
      id <- c(id, paste(ds, dis, ct, "F", sep = "/"), paste(ds, dis, ct, "M", sep = "/"))
    }
    
  }
}

all_ds_ct <- data.frame("id"=id, "count"=sex_count)
all_ds_ct <- separate(all_ds_ct, id, into = c("proj","disease", "ct", "sex"), sep = "/")

all_ds_ct$proj_dis <- paste(all_ds_ct$proj, all_ds_ct$disease, sep = "_")

write.csv(all_ds_ct, paste0(out_path, "num_cells_per_ct.csv"))

all_ds_ct <- read.csv(paste0(out_path, "num_cells_per_ct.csv"))

order_proj_dis <- c(
  "Eze_Nowakowski_2nd_trimester_Healthy",  "Velmeshev_2022_2nd trimester_Healthy", 
  "Velmeshev_2022_3rd trimester_Healthy",  "Velmeshev_2022_0-1 years_Healthy",     
  "Velmeshev_2022_1-2 years_Healthy",      "Velmeshev_2022_2-4 years_Healthy",     
  "Velmeshev_2022_10-20 years_Healthy",    "Velmeshev_2022_Adult_Healthy",         
  "DISCO_GSE157827_Healthy",  "DISCO_GSE174367_Healthy",
  "DISCO_PRJNA544731_Healthy", "DISCO_GSE157827_Alzheimer's disease",             
  "DISCO_GSE174367_Alzheimer's disease", "DISCO_PRJNA544731_Multiple Sclerosis" 
)

all_ds_ct$proj_dis <- factor(all_ds_ct$proj_dis, order_proj_dis[-1])
all_ds_ct <- all_ds_ct[order(all_ds_ct$proj_dis), ]

png(paste0(plot_path, "num_cells_per_ct.png"), height = 16, width = 10, units="in", res = 1200)
print(
  ggplot(all_ds_ct, aes(proj_dis, count, fill=proj_dis)) +
    geom_bar(stat = "identity", show.legend = T, color="black") +
    geom_hline(yintercept = 100, linetype="dashed") +
    labs( y="Number of cells", x="Datasets", fill="") +
    facet_grid(ct ~ sex, scales = "free_y") +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, color = "black"), 
      panel.spacing.x = unit(0.5, "lines"),
      strip.text.x = element_text(size=12, face="bold", colour = "black"),
      strip.text.y = element_text(size=12, face="bold", colour = "black", angle = 0),
      plot.title = element_text(size=12, face="bold", colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.title.y = element_text(size=12, face="bold", colour = "black"),
      axis.text.y = element_text(size=10, colour = "black", vjust = 0.7, hjust=0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=10, colour = "black", angle = 90),
      axis.ticks.x = element_blank(), 
      legend.position = "none"
      )
)
dev.off()
