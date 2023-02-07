library(dplyr)
library(stringr)

out_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/thesis_draft/suppl_files/"

disco <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/DISCOv1.0/DEGs_common/num_proj_sex_ct.csv")

velm <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Velmeshev_num_sex_ct_per_age.csv")
eze_nowa <- read.csv("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/outputs/Eze_Nowakowski_num_cells.csv")

disco$X <- NULL
velm$X <- NULL
eze_nowa$X <- NULL

colnames(eze_nowa) <- c("proj", "sex", "ct", "count")
colnames(velm) <- c("proj", "sex", "ct", "count")
eze_nowa$disease <- rep("Normal", nrow(eze_nowa))
velm$disease <- rep("Normal", nrow(velm))

eze_nowa <- eze_nowa %>% relocate(disease, .after = sex)
velm <- velm %>% relocate(disease, .after = sex)

velm$proj <- paste("Velmeshev_2022", velm$proj, sep = "_")
disco$proj <- paste("DISCO", disco$proj, sep = "_")


all_ds_ct <- rbind(eze_nowa, velm, disco)
all_ds_ct$sex <- str_replace_all(all_ds_ct$sex, c("Female"="F", "Male"="M"))

sex_count <- vector()
id <- vector()
for (ds in unique(all_ds_ct$proj)) {
  sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds &  all_ds_ct$sex=="F"), "count"]))
  sex_count <- c(sex_count, sum(all_ds_ct[which(all_ds_ct$proj==ds &  all_ds_ct$sex=="M"), "count"]))
  id <- c(id, paste(ds, "F", sep = "/"), paste(ds, "M", sep = "/"))
}

all_ds <- data.frame("id"=id, "count"=sex_count)
all_ds <- separate(all_ds, id, into = c("proj", "sex"), sep = "/")

write.csv(all_ds_ct, paste0(out_path, "num_cells_per_ct.csv"))
write.csv(all_ds, paste0(out_path, "num_cells.csv"))