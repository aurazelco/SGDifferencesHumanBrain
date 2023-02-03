# Author: Aura Zelco
# Brief description:
  # This script is used to find the TFs overlapping with the Cistrome DB - copy-pasted listed from the Cistrome DB on 2023/01/11 http://cistrome.org/db/#/ 
# Brief procedure:
  # 1. Reads all DEG CSV files from the SCENIC result folder
  # 2. Intersects the TFs of each files with the reference Cistrome DB
  # 3. Saves the results as a single dataframe to the output folder
# OBS: since there is a need for manual input, it is recommended to run this script in a R environment/IDE (e.g. RStudio)

#---------------------------------------------------------------------------------------------------

# sources the script containing all functions run here
source("~/Desktop/Lund_MSc/Thesis/scripts/Extra_scripts/Cistrome_DB/cistrome_DB_TFs_func.R")

# Imports the Cistrome BD TFs list and stores it as a vector
cistrome_tf <- read.table("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/Cistrome/all_TF_in_db.txt")
colnames(cistrome_tf) <- "TFs"
cistrome_tf <- subset(cistrome_tf, TFs %!in% c("All", "ATAC-seq")) # not real TFs
cistrome_tf <- as.vector(cistrome_tf$TFs)

# sets the directories where to find the GRN tsv files
main_DISCO <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/SCENIC/outputs/"
main_UCSC <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/SCENIC/"

#### UCSC
UCSC <- CistromeDB(main_UCSC, "1_GRN", F, cistrome_tf)

UCSC <- do.call(rbind, UCSC)
UCSC$sex <- str_replace_all(UCSC$sex, c("Female"="F", "Male"="M"))
UCSC <- cbind("group"=gsub("\\..*","",rownames(UCSC)),
              UCSC)
rownames(UCSC) <- NULL

for (group_id in unique(UCSC$group)) {
  print(group_id)
  for (run_id in c(1,2,3)) {
    print(setdiff(
      UCSC[which(UCSC$group==group_id & UCSC$run==run_id & UCSC$sex=="F"), "TFs_in_DB"], 
      UCSC[which(UCSC$group==group_id & UCSC$sex=="M"), "TFs_in_DB"]
    ))
    print(setdiff(
      UCSC[which(UCSC$group==group_id & UCSC$run==run_id & UCSC$sex=="M"), "TFs_in_DB"], 
      UCSC[which(UCSC$group==group_id & UCSC$sex=="F"), "TFs_in_DB"]
    ))
  }
}

#### DISCO
DISCO <- CistromeDB(main_DISCO, "1_GRN", T, cistrome_tf)

DISCO <- do.call(rbind, DISCO)
DISCO <- cbind("group"=gsub("\\..*","",rownames(DISCO)),
              DISCO)
rownames(DISCO) <- NULL

for (group_id in unique(DISCO$group)) {
  print(group_id)
  for (proj_id in unique(DISCO[which(DISCO$group==group_id), "proj"]))
  for (run_id in c(1,2,3)) {
    print(setdiff(
      DISCO[which(DISCO$group==group_id & DISCO$run==run_id & DISCO$proj==proj_id & DISCO$sex=="F"), "TFs_in_DB"], 
      DISCO[which(DISCO$group==group_id & DISCO$proj==proj_id & DISCO$sex=="M"), "TFs_in_DB"]
    ))
    print(setdiff(
      DISCO[which(DISCO$group==group_id & DISCO$run==run_id & DISCO$proj==proj_id & DISCO$sex=="M"), "TFs_in_DB"], 
      DISCO[which(DISCO$group==group_id & DISCO$proj==proj_id & DISCO$sex=="F"), "TFs_in_DB"]
    ))
  }
}

# no different TFs in any group, UCSC or DISCO