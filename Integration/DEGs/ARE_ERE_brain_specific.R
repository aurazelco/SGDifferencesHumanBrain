library(readxl)
library(ggplot2)
library(stringr)

main_int_path <- "data/Integration/"

source("scripts/Integration/Functional_analysis/Comparison_wth_ref_datasets_func.R")

tab_names <- excel_sheets(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"))
brain_tabs <- c(tab_names[which(grepl("^BRN", tab_names))], "PTTARY")

oliva <- lapply(brain_tabs, function(x) read_excel(path = paste0(main_int_path, "Oliva_2020-table-s2.xlsx"), sheet = x))
names(oliva) <- brain_tabs
oliva <- do.call(rbind, oliva)
oliva$region <- gsub("\\..*", "", rownames(oliva))
rownames(oliva) <- NULL
oliva <- as.data.frame(oliva)
oliva$chr_simplified <- str_replace_all(oliva$chr, c("chrX"="X", "chrY"="Y", "chr\\d+"="Autosome"))
oliva <- oliva[which(!is.na(oliva$HUGO_gene_id)), ]

oliva_num <- c("X"=length(unique(oliva[which(oliva$chr_simplified=="X"), "HUGO_gene_id"])),
               "Autosome" = length(unique(oliva[which(oliva$chr_simplified=="Autosome"), "HUGO_gene_id"])))
tot_genes <- 20000

oliva_num_reg <- c("BRNAMY" = length(oliva[which(oliva$region=="BRNAMY"), "HUGO_gene_id"]), 
                   "BRNACC"  = length(oliva[which(oliva$region=="BRNACC"), "HUGO_gene_id"]),
                   "BRNCDT"  = length(oliva[which(oliva$region=="BRNCDT"), "HUGO_gene_id"]),
                   "BRNCHB"  = length(oliva[which(oliva$region=="BRNCHB"), "HUGO_gene_id"]),
                   "BRNCHA"  = length(oliva[which(oliva$region=="BRNCHA"), "HUGO_gene_id"]),
                   "BRNCTXA" = length(oliva[which(oliva$region=="BRNCTXA"), "HUGO_gene_id"]),
                   "BRNCTXB" = length(oliva[which(oliva$region=="BRNCTXB"), "HUGO_gene_id"]),
                   "BRNHPP"  = length(oliva[which(oliva$region=="BRNHPP"), "HUGO_gene_id"]),
                   "BRNHPT"  = length(oliva[which(oliva$region=="BRNHPT"), "HUGO_gene_id"]),
                   "BRNNCC"  = length(oliva[which(oliva$region== "BRNNCC"), "HUGO_gene_id"]),
                   "BRNPTM"  = length(oliva[which(oliva$region=="BRNPTM"), "HUGO_gene_id"]),
                   "BRNSPC"  = length(oliva[which(oliva$region=="BRNSPC"), "HUGO_gene_id"]),
                   "BRNSNG" = length(oliva[which(oliva$region=="BRNSNG"), "HUGO_gene_id"]),
                   "PTTARY" = length(oliva[which(oliva$region=="PTTARY"), "HUGO_gene_id"])
)


# On all regions together

main_DISCO <- "DISCOv1.0/DEGs_proj_adjust_pval/"
main_UCSC <- "UCSC/DEGs_adjust_pval/"

# Vectors to save the different sub-groups of DISCO and UCSC
# the first folder "exta_files" is excluded
sub_projs <- list.dirs(main_DISCO, full.names = F, recursive = F)[-1]
sub_UCSC <- list.dirs(main_UCSC, full.names = F, recursive = F)[-1]

# Import all the CSVs from the different ages/groups - slightly different file tree structure requires a different approach for UCSC

disco <- ImportDatasets(main_DISCO, sub_projs, UCSC_flag = "no", individual_projs = T)
names(disco[[1]]) <- str_replace_all(names(disco[[1]]), "Normal", "Healthy")
UCSC <- ImportDatasets(main_UCSC, sub_UCSC, UCSC_flag = "yes", individual_projs = F)

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
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
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



# Combines all dataframes into one df
sexes <- CreateSexDf(c(disco[[1]], UCSC[[1]]), unified_annotation)
groups_order <- c(
  "Velmeshev_2022_2nd_trimester",           
  "Velmeshev_2022_3rd_trimester", 
  "Velmeshev_2022_0_1_years",                
  "Velmeshev_2022_1_2_years",            
  "Velmeshev_2022_2_4_years",  
  "Velmeshev_2022_10_20_years",      
  "Velmeshev_2022_Adult",
  "Healthy_GSE157827",              
  "Healthy_GSE174367",               
  "Healthy_PRJNA544731", 
  "Alzheimer's disease_GSE157827",
  "Alzheimer's disease_GSE174367",
  "Multiple Sclerosis_PRJNA544731" 
)
count_oliva_all <- CountOliva(main_int_path, sexes, oliva, groups_order)
count_oliva_all$oliva_perc <- count_oliva_all$oliva_count * 100 / count_oliva_all$tot_degs_count


# ARE and ERE refs

ARE <- read_excel("data/DISCOv1.0/DEGs_proj_adjust_pval/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("data/DISCOv1.0/DEGs_proj_adjust_pval/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

ARE_full <- c()
reg <- c()
chr <- c()
for (region in unique(oliva$region)) {
  for (chrom in unique(oliva[which(oliva$region==region), "chr_simplified"])) {
    ARE_full <- c(ARE_full, sum(unique(ARE$fullsites) %in% oliva[which(oliva$region==region & oliva$chr_simplified==chrom), "HUGO_gene_id"]))
    reg <- c(reg, region)
    chr <- c(chr, chrom)
}
}

ARE_half <- c()
for (region in unique(oliva$region)) {
  for (chrom in unique(oliva[which(oliva$region==region), "chr_simplified"])) {
    ARE_half <- c(ARE_half, sum(unique(ARE$halfsites) %in% oliva[which(oliva$region==region & oliva$chr_simplified==chrom), "HUGO_gene_id"]))
    reg <- c(reg, region)
    chr <- c(chr, chrom)
  }
}


ERE_full <- c()
for (region in unique(oliva$region)) {
  for (chrom in unique(oliva[which(oliva$region==region), "chr_simplified"])) {
    ERE_full <- c(ERE_full, sum(unique(ERE$`Hs Gene Name`) %in% oliva[which(oliva$region==region & oliva$chr_simplified==chrom), "HUGO_gene_id"]))
    reg <- c(reg, region)
    chr <- c(chr, chrom)
  }
}

sites <- c(rep("ARE_full", length(ARE_full)), rep("ARE_half", length(ARE_half)), rep("ERE", length(ERE_full)))
count_re <- c(ARE_full, ARE_half, ERE_full)

re <- data.frame(sites, count_re, reg, chr)


ggplot(re, aes(x= sites, y= count_re, fill=reg)) +
  geom_bar(position="dodge", stat = "identity")



ARE_brain <- sum(unique(c(ARE$fullsites, ARE$halfsites)) %in% oliva$HUGO_gene_id) * 100 / length(unique(c(ARE$fullsites, ARE$halfsites)))
ERE_brain <- sum(unique(ERE$`Hs Gene Name`) %in% oliva$HUGO_gene_id) * 100 / length(unique(ERE$`Hs Gene Name`))

#ARE_brain
#[1] 5.347537
#ERE_brain
#[1] 2.221334

library(ggVennDiagram) 
library(ggpubr)

# https://stackoverflow.com/questions/70228591/set-the-color-of-categories-in-venn-diagram-in-r
plot_venn <- function (x, show_intersect, set_color, set_size, label, label_geom, 
                       label_alpha, label_color, label_size, label_percent_digit, 
                       label_txtWidth, edge_lty, edge_size, ...)  {
  venn <- Venn(x)
  data <- process_data(venn)
  p <- ggplot() + geom_sf(aes_string(fill = "count"), data = data@region) + 
    geom_sf(aes_string(color = "name"), data = data@setEdge, 
            show.legend = F, lty = edge_lty, size = edge_size, color = set_color) + 
    geom_sf_text(aes_string(label = "name"), data = data@setLabel, 
                 size = set_size, color = set_color) + theme_void()
  if (label != "none" & show_intersect == FALSE) {
    region_label <- data@region %>% dplyr::filter(.data$component == 
                                                    "region") %>% dplyr::mutate(percent = paste(round(.data$count * 
                                                                                                        100/sum(.data$count), digits = label_percent_digit), 
                                                                                                "%", sep = "")) %>% dplyr::mutate(both = paste(.data$count, 
                                                                                                                                               paste0("(", .data$percent, ")"), sep = "\n"))
    if (label_geom == "label") {
      p <- p + geom_sf_label(aes_string(label = label), 
                             data = region_label, alpha = label_alpha, color = label_color, 
                             size = label_size, lineheight = 0.85, label.size = NA)
    }
    if (label_geom == "text") {
      p <- p + geom_sf_text(aes_string(label = label), 
                            data = region_label, alpha = label_alpha, color = label_color, 
                            size = label_size, lineheight = 0.85)
    }
  }
  if (show_intersect == TRUE) {
    items <- data@region %>% dplyr::rowwise() %>% dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, 
                                                                                                collapse = " "), width = label_txtWidth)) %>% sf::st_as_sf()
    label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
    p <- ggplot(items) + geom_sf(aes_string(fill = "count")) + 
      geom_sf_text(aes_string(label = "name"), data = data@setLabel, 
                   inherit.aes = F) + geom_text(aes_string(label = "count", 
                                                           text = "text"), x = label_coord[, 1], y = label_coord[, 
                                                                                                                 2], show.legend = FALSE) + theme_void()
    ax <- list(showline = FALSE)
    p <- plotly::ggplotly(p, tooltip = c("text")) %>% plotly::layout(xaxis = ax, 
                                                                     yaxis = ax)
  }
  p
}
assignInNamespace(x="plot_venn", value=plot_venn, ns="ggVennDiagram")  

out_path <- paste0(main_int_path, "ARE_ERE_brain/")
dir.create(out_path)


all_plt <-  ggVennDiagram(x= list(unique(ARE$fullsites, ARE$halfsites), ERE$`Hs Gene Name`, unique(oliva$HUGO_gene_id)),
                category.names = c("ARE", "ERE", "Brain"),
                label_alpha = 100, 
                color =  c("ARE" = "#F8766D", "ERE" ="#00BA38",'Brain' = "#619CFF"),
                set_color = c("ARE" = "#F8766D", "ERE" ="#00BA38",'Brain' = "#619CFF")) +
    theme(legend.position = "none")



ARE_plt <- ggVennDiagram(x= list(unique(ARE$fullsites, ARE$halfsites), unique(oliva$HUGO_gene_id)),
                category.names = c("ARE", "Brain"),
                label_alpha = 100, 
                label_percent_digit = 2,
                color =  c("ARE" = "#F8766D", 'Brain' = "#619CFF"),
                set_color = c("ARE" = "#F8766D", 'Brain' = "#619CFF")) +
    theme(legend.position = "none")


ERE_plt <- ggVennDiagram(x= list(ERE$`Hs Gene Name`, unique(oliva$HUGO_gene_id)),
                category.names = c("ERE", "Brain"),
                label_alpha = 100, 
                label_percent_digit = 2,
                color =  c("ERE" ="#00BA38",'Brain' = "#619CFF"),
                set_color = c("ERE" ="#00BA38",'Brain' = "#619CFF")) +
    theme(legend.position = "none")


pdf(paste0(out_path, "Venn_diagram_combined.pdf"), width = 8)
print(ggarrange(ARE_plt, ERE_plt, 
                labels = c("A", "B"),
                nrow=2,
                align = "v"))
dev.off()