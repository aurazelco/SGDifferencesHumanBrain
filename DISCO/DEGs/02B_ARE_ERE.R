main <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/outputs/"

source("/Users/aurazelco/Desktop/Lund_MSc/Thesis/scripts/DISCO/DEGs/02B_ARE_ERE_func.R")

####### MAIN

sub_disease <- list.dirs(main, recursive=FALSE, full.names = FALSE)

ARE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/AREsitesHuman.xlsx",
                  skip=1)
colnames(ARE) <- c("fullsites", "halfsites")

ERE <- read_excel("Desktop/Lund_MSc/Thesis/data/DISCOv1.0/20220817_DEGs/extra_files/Hs_allEREs.xls")
EREgene <- ERE$`Hs Gene Name`

# NORMAL
normal <- sub_disease[3]
AnalysisARE_ERE(main, normal, ARE, EREgene)

#AD
ad <- sub_disease[1]
AnalysisARE_ERE(main, ad, ARE, EREgene)
