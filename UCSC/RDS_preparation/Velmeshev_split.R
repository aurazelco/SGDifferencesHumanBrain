library(Seurat)
library(SeuratObject)
library(data.table)
library(stringr)
library(stringi)
library(tidyr)
library(ggplot2)
library(matrixStats)
library(dplyr)
`%!in%` <- Negate(`%in%`)

rds_path <- "UCSC/Seurat_UCSC/"
main_deg <- "UCSC/DEGs/"
main_scenic <- "UCSC/SCENIC/"

main <- "UCSC/outputs/Velmeshev/"

# Clusters annotation using Fig 1B and Fig S1a from paper

ann_clusters <- list("Astrocytes" = c(3),
                  "Microglia" = c(25),
                  "OPCs" = c(9),
                  "Oligodendrocytes" = c(11),
                  "Excitatory neurons" = c(0,13,4,7,21,12,19,24,6,28,23,10,22,26,15, 30),
                  "Dorsal progenitors" = c(18),
                  "Ventral progenitors" = c(5),
                  "Interneurons" = c(17, 31,32,8,16, 14),
                  "Vascular cells" = c(20),
                  "Debris" = c(1,2),
                  "Unknown" = c(29, 27)
                  )

Reduce(intersect, ann_clusters)

cts <- vector()
og_clusters <- vector()

for (ct in names(ann_clusters)) {
  cts <- c(cts, rep(ct, length(ann_clusters[[ct]])))
  og_clusters <- c(og_clusters, ann_clusters[[ct]])
}
ann_df <- data.frame(cts, og_clusters)


# from this tutorial UCSC CellBrowser https://cellbrowser.readthedocs.io/en/master/load.html

####################################################################################################
#
# VELMESHEV 2022
#
####################################################################################################

####  Separate cell barcodes by age group - prepare indexes for file splitting

meta <- read.csv("UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
split_barcodes <- list()
for (age in unique(meta$age)) {
  barcodes_id <- meta[which(meta$age==age), "cell"]
  #split_barcodes <- append(split_barcodes, list((paste(c(age, "gene", barco(des_id), collapse = ","))))
  split_barcodes <- append(split_barcodes, list((paste(c("gene", barcodes_id), collapse = ","))))
}
Reduce(intersect, split_barcodes)

sink("UCSC/UCSC_downloads/Velmeshev_barcodes.txt")
print(split_barcodes)
sink()


names(split_barcodes) <- str_replace_all(unique(meta$age), c(" "="_", "-"="_"))


meta <- read.csv("UCSC/UCSC_downloads/new_meta_Velmeshev_2022.csv", header=T, sep=",", as.is=T, row.names=1)
split_index <- list()
index_names <- c()
problematic <- c("0-1 years", "3rd trimester", "2nd trimester")
for (age in unique(meta$age)) {
  indexes <- sort(as.numeric(rownames(meta[which(meta$age==age),])))
  if (age %in% problematic) {
    indexes <- indexes + 1
    half <- length(indexes)/2
    if (half %% 1 == 0) {
      ind1 <- indexes[1:half]
      ind2 <- indexes[(half+1):(length(indexes))]
    } else {
      ind1 <- indexes[1:(half - 0.5)]
      ind2 <- indexes[(half + 0.5):(length(indexes))]
    }
    split_index <- append(split_index, list(c(1,ind1)))
    split_index <- append(split_index, list(c(1,ind2)))
    ind1_name <- paste0(str_replace_all(age, c(" "="_", "-"="_")), "_1")
    index_names <- c(index_names, ind1_name)
    ind2_name <- paste0(str_replace_all(age, c(" "="_", "-"="_")), "_2")
    index_names <- c(index_names, ind2_name)
  } else {
    #split_index <- append(split_index, list(paste((rownames(meta[which(meta$age==age),])), collapse = ",")))
    indexes <- indexes + 1
    split_index <- append(split_index, list(c(1,indexes)))
    index_names <- c(index_names, str_replace_all(age, c(" "="_", "-"="_")))
  }
}
names(split_index) <- index_names
rm(age, half, ind1, ind1_name, ind2, ind2_name, index_names, indexes)

Reduce(intersect, split_index)

lapply(1:length(names(split_index)), function(x) write.table(
  t(split_index[[x]]),
  paste0("UCSC/UCSC_downloads/Velmeshev_split/", names(split_index)[x], ".txt"),
  sep=",",
  eol="",
  row.names = F,
  col.names = F))

probs_files <- list.files("UCSC/UCSC_downloads/Velmeshev_problematic", pattern = "*.txt", full.names = T)
probs_names <- list.files("UCSC/UCSC_downloads/Velmeshev_problematic", pattern = "*.txt", full.names = F)
probs_names <- str_remove_all(probs_names, ".txt")
mat_filt <- fread("UCSC/UCSC_downloads/exprMtx_filt_Velmeshev_2022.tsv.gz")

for (i in 1:length(probs_files)) {
  prob <- read.table(probs_files[i], sep=",", header = F)
  prob <- as.numeric(prob)
  mat_prob <- mat_filt[,..prob]
  gz_prob <- gzfile(paste0("UCSC/UCSC_downloads/Velmeshev_outs/", probs_names[i], ".tsv.gz"), "w")
  write.table(mat_prob, gz_prob, sep="\t")
  close(gz_prob)
}


