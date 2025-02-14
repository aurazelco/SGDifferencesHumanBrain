library(readxl)
library(ggplot2)
library(stringr)
library(tidyverse)

main_int_path <- "data/Integration/"

source("scripts/Integration/Functional_analysis/Comparison_wth_ref_datasets_func.R")

tab_names <- excel_sheets(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"))
tab_names <- tab_names[-c(1,2)]

oliva <- lapply(tab_names, function(x) read_excel(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"), sheet = x))
names(oliva) <- tab_names
oliva <- do.call(rbind, oliva)
oliva$tissue <- gsub("\\..*", "", rownames(oliva))
rownames(oliva) <- NULL

oliva <- oliva %>%
          mutate(gene_symbol=ifelse(is.na(HUGO_gene_id), ENSEMBL_gene_id, HUGO_gene_id)) %>%
          select(tissue, gene_symbol)

brain_regions <- c(tab_names[which(grepl("^BRN", tab_names))], "PTTARY")

oliva <- oliva %>% 
          bind_rows(oliva %>%
                      filter(tissue %in% brain_regions) %>%
                      mutate(tissue="BRN_TOTAL") %>%
                      unique()) %>%
          arrange(tissue)


# ARE and ERE refs

ARE_genes <- read_excel("DISCOv1.0/DEGs_proj_adjust_pval/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE_genes) <- c("fullsites", "halfsites")

ARE_genes <- unique(ARE_genes$fullsites, ARE_genes$halfsites)

ERE_genes <- read_excel("DISCOv1.0/DEGs_proj_adjust_pval/extra_files/Hs_allEREs.xls")
ERE_genes <- unique(ERE_genes$`Hs Gene Name`)

oliva <- oliva %>%
                mutate(ARE_common = ifelse(gene_symbol %in% ARE_genes, 1, 0),
                       ERE_common = ifelse(gene_symbol %in% ERE_genes, 1, 0))

ARE_common <- oliva %>%
              filter(ARE_common==1) %>%
              count(tissue, ARE_common) %>%
              rename(ARE_found=n) %>%
              select(tissue, ARE_found)

ERE_common <- oliva %>%
              filter(ERE_common==1) %>%
              count(tissue, ERE_common) %>%
              rename(ERE_found=n) %>%
              select(tissue, ERE_found)


oliva_perc <- Reduce(function(...) merge(..., by='tissue', all.x=TRUE), list(oliva %>%
                                                                                count(tissue) %>%
                                                                                rename(total=n),
                                                                              ARE_common,
                                                                              ERE_common)) %>%
  mutate(ARE = ARE_found * 100 / total,
         ERE =ERE_found * 100 / total) %>%
  pivot_longer(cols=c(ARE, ERE),
               names_to = "RE",
               values_to ="percent") %>%
  select(-c(total, ARE_found, ERE_found)) 



pdf("Integration/ARE_ERE_organs/Faceted_ARE_ERE_organ_specificity.pdf")
print(ggplot(oliva_perc, aes(tissue, percent, fill=RE)) +
        geom_bar(stat="identity", color="black") +
        labs(x="Tissues", y="% of RE sites", fill="Response element sites") +
        ylim(0, 100) +
        facet_wrap(~RE, nrow = 2) +
        scale_fill_manual(values = c('ARE'="#E1AD01", 'ERE'="#39B600")) +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              panel.background = element_blank(), 
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(size=12, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
              axis.title.x=element_text(size=14, face="bold", colour = "black"),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title.y = element_text(size=14, face="bold", colour = "black"),
              legend.position = "bottom", 
              legend.title = element_text(size=14, face="bold", colour = "black"),
              legend.text = element_text(size=14, face="bold", colour = "black"),
              strip.text.x = element_text(size=14, face="bold", colour = "black"),
              strip.text.y.right = element_text(size=14, face="bold", colour = "black", angle = 0))
)
dev.off()


