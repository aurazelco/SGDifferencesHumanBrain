library(ggplot2)
library(ggpubr)
`%!in%` <- Negate(`%in%`)


pubmed_all <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Thesis_draft/PubMed_Timeline_Results_by_Year.csv", sep=",", row.names = NULL)
pubmed_all <- pubmed_all[-c(1,2),]
colnames(pubmed_all) <- c("Year", "Count_all")
str(pubmed_all)
pubmed_all$Year <- as.numeric(pubmed_all$Year)
pubmed_all$Count_all <- as.numeric(pubmed_all$Count_all)
pubmed_all <- pubmed_all[order(pubmed_all$Year),]

pubmed_dev <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/Thesis_draft/PubMed_Timeline_Results_by_Year_development.csv", sep=",", row.names = NULL)
pubmed_dev <- pubmed_dev[-c(1,2),]
colnames(pubmed_dev) <- c("Year", "Count_dev")
str(pubmed_dev)
pubmed_dev$Year <- as.numeric(pubmed_dev$Year)
pubmed_dev$Count_dev <- as.numeric(pubmed_dev$Count_dev)
pubmed_dev <- pubmed_dev[order(pubmed_dev$Year),]

missing_years <- pubmed_all[which(pubmed_all$Year %!in% pubmed_dev$Year), "Year"]

missing_df <- data.frame(missing_years, rep(0, length(missing_years)))
colnames(missing_df) <- c("Year", "Count_dev")

pubmed_dev <- rbind(missing_df, pubmed_dev)
pubmed_dev <- pubmed_dev[order(pubmed_dev$Year),]


pubmed <- merge(pubmed_all, pubmed_dev, by="Year", all.x=TRUE, all.y=TRUE)
pubmed$adult <- pubmed$Count_all - pubmed$Count_dev
pubmed <- reshape2::melt(pubmed[c(1, 3:4)], "Year")
pubmed$Year <- as.factor(pubmed$Year)
colnames(pubmed) <- c("year", "age", "count")
levels(pubmed$age) <- c("development", "non-development")

all <- ggplot(pubmed, aes(year, count, fill=age)) +
  geom_bar(stat="identity", color="black") +
  scale_fill_manual(values = c("non-development"= "#D55E00", "development" = "#009E73")) +
  labs(x="Year", y="Publications", fill="Age", title="A)") +
  scale_x_discrete(breaks=seq(min(as.numeric(levels(pubmed$year))),
                              max(as.numeric(levels(pubmed$year))),
                              10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.title = element_text(size=12, face="bold", colour = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold", colour = "black"),
        legend.position = "bottom", 
        legend.title = element_text(size=14, face="bold", colour = "black"))

perc <- ggplot(pubmed, aes(year, count, fill=age)) +
  geom_bar(stat="identity", color="black", position = "fill") +
  scale_fill_manual(values = c("non-development"= "#D55E00", "development" = "#009E73")) +
  labs(x="Year", y="Publications (%)", fill="Age", title="B)") +
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
ggsave("/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/images/publications.pdf")

