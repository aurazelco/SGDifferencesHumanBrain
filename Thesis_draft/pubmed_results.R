# Author: Aura Zelco
# Brief description:
  # This script is used to generate bar plots based on PubMed publication data
# Brief procedure:
  # 1. Reads the corresponding csv files (in this case 2)
  # 2. merges the dataframes so both information are present in the final csv
  # 3. generate the plots
  # 4. generate a composite plot
  # 5. saves said pot as PDF

# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. Rstudio)


library(ggplot2)
library(ggpubr)
`%!in%` <- Negate(`%in%`)

files_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Thesis_draft/extra_files/"
dir.create("/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/images/", showWarnings = F, recursive = T)

sexdiff_brain_all <- read.csv(paste0(files_path, "PubMed_Timeline_Results_by_Year_sexdiff_AND_brain.csv"), 
                              sep=",", row.names = NULL, skip = 1)
colnames(sexdiff_brain_all) <- c("Year", "Count_all")
sexdiff_brain_all <- sexdiff_brain_all[order(sexdiff_brain_all$Year),]

sexdiff_brain_dev <- read.csv(paste0(files_path, "PubMed_Timeline_Results_by_Year_sexdiff_AND_brain_AND_dev.csv"), 
                       sep=",", row.names = NULL, skip = 1)
colnames(sexdiff_brain_dev) <- c("Year", "Count_dev")
sexdiff_brain_dev <- sexdiff_brain_dev[order(sexdiff_brain_dev$Year),]

missing_years <- sexdiff_brain_all[which(sexdiff_brain_all$Year %!in% sexdiff_brain_dev$Year), "Year"]

missing_df <- data.frame(missing_years, rep(0, length(missing_years)))
colnames(missing_df) <- c("Year", "Count_dev")

sexdiff_brain_dev <- rbind(missing_df, sexdiff_brain_dev)
sexdiff_brain_dev <- sexdiff_brain_dev[order(sexdiff_brain_dev$Year),]


pubmed <- merge(sexdiff_brain_all, sexdiff_brain_dev, by="Year", all.x=TRUE, all.y=TRUE)
pubmed$non_dev <- pubmed$Count_all - pubmed$Count_dev
pubmed <- reshape2::melt(pubmed[c(1, 3:4)], "Year")
pubmed$Year <- as.factor(pubmed$Year)
colnames(pubmed) <- c("year", "publ", "count")
levels(pubmed$publ) <- c("Development", "Non-development")

all <- ggplot(pubmed, aes(year, count, fill=publ)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("Non-development"= "#D55E00", "Development" = "#009E73")) +
  labs(x="Year", y="Publications", fill="Publications type", title="A)") +
  scale_x_discrete(breaks=seq(min(as.numeric(levels(pubmed$year))),
                              max(as.numeric(levels(pubmed$year))),
                              10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.x = element_text(size=12, colour = "black", vjust = 0.7, hjust=0.5),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=14, face="bold", colour = "black"))

perc <- ggplot(pubmed, aes(year, count, fill=publ)) +
  geom_bar(stat="identity", color="black", position = "fill") +
  scale_fill_manual(values = c("Non-development"= "#D55E00", "Development" = "#009E73")) +
  labs(x="Year", y="Publications (%)", fill="Publications type", title="B)") +
  scale_x_discrete(breaks=seq(min(as.numeric(levels(pubmed$year))),
                              max(as.numeric(levels(pubmed$year))),
                              10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=14, face="bold", colour = "black"),
        axis.text.x = element_text(size=12, colour = "black", vjust = 0.7, hjust=0.5),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=14, face="bold", colour = "black"))


ggarrange(all, perc, nrow = 2, common.legend = T, legend = "bottom")
ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/images/publications.png")

