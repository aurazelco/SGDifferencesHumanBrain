library(ggplot2)
library(ggpubr)
library(stringr)
#library(scales)

main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/"


####################################################################################################
#
# PLOTS
#
####################################################################################################


########### NUM OF SAMPLES AND CELLS PER SAMPLE, DIVIDED BY AGE AND SEX

num_cells <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_cells_per_age.csv")
num_cells[,1] <- NULL

trim_2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_2022_2nd_trimester_num_cells.csv")
trim_2[,1] <- NULL
yo_10_20 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_2022_10_20_years_num_cells.csv")
yo_10_20[,1] <- NULL

num_cells <- rbind(num_cells, trim_2, yo_10_20)

rm(trim_2, yo_10_20)

num_cells$age <- factor(num_cells$age, c(
  "2nd trimester", "3rd trimester","0-1 years","1-2 years","2-4 years","4-10 years","10-20 years", "Adult")
)

# needed it to correct it only once
#num_cells[which(num_cells$id=="1-1547-BA24"), "sex"] <- "Female"

p1 <- ggplot(num_cells, aes(age, fill=sex)) +
  geom_bar(position = "dodge") +
  scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  labs(x="Age", y="Number of samples", fill="Sex") +
  #facet_wrap(~proj, scales = "free", nrow=1, drop = TRUE) +
  geom_hline(yintercept = 3, linetype=2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

p2 <- ggplot(num_cells, aes(id, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  #scale_y_continuous(breaks= seq(0, nrow(num_cells),by=1)) +
  geom_hline(yintercept = 500, linetype = "dashed") +
  labs(x="Age", y="Number of cells/sample", fill="Sex") +
  facet_wrap(~age, scales = "free", nrow=2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

num_cells_proj <- ggarrange(p1, p2, common.legend = T, legend = "bottom", nrow = 2)

pdf(paste0(main, "Velmeshev_num_samples_and_cells_per_age.pdf"), height = 15)
print(num_cells_proj)
dev.off()

write.csv(num_cells, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_cells_per_age.csv")


########### NUM OF CELLS PER CLUSTER AND SEX, DIVIDED BY AGE


num_sex_ct <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")
num_sex_ct[,1] <- NULL

trim_2 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_2022_2nd_trimester_num_sex_ct_per_age.csv")
trim_2[,1] <- NULL
yo_10_20 <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_2022_10_20_years_num_sex_ct_per_age.csv")
yo_10_20[,1] <- NULL

num_sex_ct <- rbind(num_sex_ct, trim_2, yo_10_20)

rm(trim_2, yo_10_20)

num_sex_ct$proj <- str_remove_all(num_sex_ct$proj, "Velmeshev_2022_")
num_sex_ct$proj <- str_replace_all(num_sex_ct$proj, c("_years" = " years",  "_trimester" = " trimester", "_" = "-"))


num_sex_ct$proj <- factor(num_sex_ct$proj, c(
  "2nd trimester", "3rd trimester","0-1 years","1-2 years","2-4 years","4-10 years","10-20 years", "Adult")
)

# needed it to correct it only once
#num_sex_ct[which(num_sex_ct$id=="1-1547-BA24"), "sex"] <- "Female"

p3 <- ggplot(num_sex_ct, aes(ct, Freq, fill=sex)) +
  geom_bar(stat="identity", position="dodge") +
  #scale_y_continuous(breaks= seq(0, nrow(num_sex_ct),by=1)) +
  geom_hline(yintercept = 500, linetype = "dashed") +
  labs(x="Cell types", y="Number of cells/cell type", fill="Sex") +
  facet_wrap(~proj, scales = "free", nrow=2) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=12, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, colour = "black",angle = 90, vjust = 0.7, hjust=0.5),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size=12, face="bold", colour = "black"),
        axis.text.y = element_text(size=10, colour = "black",angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = "bottom", 
        legend.title = element_text(size=12, face="bold", colour = "black"),
        strip.text = element_text(size=12, face="bold", colour = "black"))

pdf(paste0(main, "Velmeshev_num_cells_per_sex_ct.pdf"))
print(p3)
dev.off()

write.csv(num_sex_ct, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev/Velmeshev_num_sex_ct_per_age.csv")


