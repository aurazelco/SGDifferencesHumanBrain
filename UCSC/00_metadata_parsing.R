library(stringr)
library(stringi)
library(tidyr)
library(pdftools)
library(readxl)


input_path <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads"

all_meta_files <- list.files(input_path, pattern = "\\.tsv$", full.names = T)

all_meta_files <- lapply(all_meta_files, function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})
names(all_meta_files) <- list.files(input_path, pattern = "\\.tsv$", full.names = F)
names(all_meta_files) <- str_remove_all(names(all_meta_files), "meta_")
names(all_meta_files) <- str_remove_all(names(all_meta_files), ".tsv")


############## vectors to summarize all info

age <- vector()
cell_type <- vector()
area <- vector()
sex  <- vector()
num_samples <- vector()

############## "Bhaduri_2021_neocortex"   

colnames(all_meta_files[[1]])

all_meta_files[[1]]$Samples <- stri_replace_last_fixed(all_meta_files[[1]]$Cell.Name, "_", "")
all_meta_files[[1]]$Samples <- stri_replace_last_fixed(all_meta_files[[1]]$Samples, "_", ".")
all_meta_files[[1]]$Samples <- gsub("\\..*", "", all_meta_files[[1]]$Samples)
all_meta_files[[1]]$Samples <- str_replace_all(all_meta_files[[1]]$Samples, c("_motor"="", "_PFC"="", "_V1"=""))
all_meta_files[[1]]$Samples <- str_replace_all(all_meta_files[[1]]$Samples, "gw", "GW")

age <- c(age, paste(unique(all_meta_files[[1]]$Age),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[1]]$ConsensusCellType...Final),collapse=", "))
area <- c(area, paste(unique(all_meta_files[[1]]$Area),collapse=", "))
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[1]]$Samples)))

############## "Bhaduri_2021_whole_brain" 

colnames(all_meta_files[[2]])

all_meta_files[[2]]$age <- paste0("GW", all_meta_files[[2]]$age)

age <- c(age, paste(unique(all_meta_files[[2]]$age),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[2]]$cell_type),collapse=", "))
area <- c(area, paste(unique(all_meta_files[[2]]$area),collapse=", "))
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[2]]$individual)))

############## "Darmanis_15"              

colnames(all_meta_files[[3]])

age <- c(age, paste(unique(all_meta_files[[3]]$GEO_Sample_age),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[3]]$biosample_cell_type),collapse=", "))
area <- c(area, paste(unique(all_meta_files[[3]]$sample_category),collapse=", "))
sex  <- c(sex, "no data")
num_samples <- c(num_samples, "no data")

############## "Eze_2021"  

colnames(all_meta_files[[4]])

all_meta_files[[4]]$Age_.Carnegie_Stage. <- paste0("GW", all_meta_files[[4]]$Age_.Carnegie_Stage.)

age <- c(age, paste(unique(all_meta_files[[4]]$Age_.Carnegie_Stage.),collapse=", "))
cell_type <- c(cell_type, "no data")
area <- c(area, paste(unique(all_meta_files[[4]]$Area_As_Annotated),collapse=", "))
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[4]]$Individual)))

############## "Han_2020"   

colnames(all_meta_files[[5]])

age <- c(age, "no data")
cell_type <- c(cell_type,  paste(unique(all_meta_files[[5]][which(all_meta_files[[5]]$sample=="FetalBrain"), "cell_type"]),collapse=", "))
area <- c(area, paste(unique(all_meta_files[[5]]$sample),collapse=", "))
sex  <- c(sex, "no data") 
num_samples <- c(num_samples, length(unique(all_meta_files[[5]][which(all_meta_files[[5]]$organ=="Brain"), "donor"])))

############## "Muller_2018"      

colnames(all_meta_files[[6]])

age <- c(age, "no data")
cell_type <- c(cell_type, paste(unique(all_meta_files[[6]][which(all_meta_files[[6]]$Tumor.Normal.Classification=="Normal"), "Cell.Type.Assignment"]),collapse=", "))
area <- c(area, "no data")
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[6]][which(all_meta_files[[6]]$Tumor.Normal.Classification=="Normal"), "Study_ID"])))


############## "Nowakowski_2017"          

colnames(all_meta_files[[7]])

all_meta_files[[7]]$Age_in_Weeks <- paste0("GW", all_meta_files[[7]]$Age_in_Weeks)

age <- c(age, paste(unique(all_meta_files[[7]]$Age_in_Weeks),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[7]]$WGCNAcluster),collapse=", "))
area <- c(area, paste(unique(all_meta_files[[7]]$RegionName),collapse=", "))
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[7]]$Name)))


############## "unknown"             

colnames(all_meta_files[[8]])

all_meta_files[[8]]$samples <- gsub(".*-","",all_meta_files[[8]]$cell)

age <- c(age, paste(unique(all_meta_files[[8]]$age),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[8]]$cell_type),collapse=", "))
area <- c(area, "no data")
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[8]][which(all_meta_files[[8]]$age %in%  c("3rd trimester", "2nd trimester", "0-1 years" )), "samples"])))

############## "van_Bruggen_2022"    

colnames(all_meta_files[[9]])

age <- c(age, paste(unique(all_meta_files[[9]]$Age),collapse=", "))
cell_type <- c(cell_type, paste(unique(all_meta_files[[9]]$Clusters),collapse=", "))
area <- c(area, "no data")
sex  <- c(sex, "no data")
num_samples <- c(num_samples, length(unique(all_meta_files[[9]]$Sample)))

############## "Velmeshev_2022"   

colnames(all_meta_files[[10]])

all_meta_files[[10]]$samples <- gsub(".*-","",all_meta_files[[10]]$cell)

age <- c(age, paste(unique(all_meta_files[[10]]$age),collapse=", "))
cell_type <- c(cell_type, "no data")
area <- c(area, "no data")
sex  <- c(sex, paste(unique(all_meta_files[[10]]$sex),collapse=", "))
num_samples <- c(num_samples, length(unique(all_meta_files[[10]][which(all_meta_files[[10]]$age %in%  c("3rd trimester", "2nd trimester", "0-1 years" )), "samples"])))

# unknown may be Velmeshev_2022

############## create df

df_meta <- data.frame(names(all_meta_files), age, num_samples, sex, area, cell_type)
colnames(df_meta) <- c("projects", "age", "num_samples", "sex", "area", "cell_type")

# time to scavenge for extra metadata

# Muller actually does not have fetal data it seems -> discuss with Anagha

# no additional metadata for Van Bruggen, but it is mentioned in the paper that the samples are from the forebrain
df_meta[which(df_meta$projects=="van_Bruggen_2022"), "area"] <- "forebrain"

# no additional metadata for Velmeshev

# put suppl together
suppl_files <- list()
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Bhaduri_2021_suppl.xlsx",
                                                    sheet = 5))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Bhaduri_2021_suppl.xlsx",
                                 sheet = 1))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Darmanis_2015_suppl.xlsx",
                                                    sheet = 1, col_names = T, skip=1))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Eze_2021_suppl.xlsx",
                                                    sheet = 6))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Han_2020_suppl.xlsx",
                                                    sheet = 2, skip = 2))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Muller_2018_suppl.xlsx"))))
suppl_files <- append(suppl_files,
                      list(as.data.frame(read_excel("/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/UCSC_downloads/meta_Nowakowski_2017_suppl.xlsx",
                                                    sheet = 2))))
names(suppl_files) <- c("meta_Bhaduri_2021_neocortex_suppl.xlsx", list.files(input_path, pattern = "\\.xlsx$"))
names(suppl_files) <- str_remove_all(names(suppl_files), "meta_")
names(suppl_files) <- str_remove_all(names(suppl_files), "_suppl.xlsx")
names(suppl_files)[2] <- paste0(names(suppl_files)[2], "_whole_brain")
colnames(suppl_files[["Han_2020"]]) <- c("Tissues/Cell lines", "Donor",                                              
                                           "Age", "Gender",                                              
                                           "Abortion/DCD(Donation after Cardiac Death)/Operation", "Cause of death",                                      
                                           "Medical history", "Sampling location",                                   
                                           "Dissociation reagent", "Reagent concentration",                                
                                           "Enzymatic digestion time" )


df_meta[which(df_meta$projects=="Darmanis_2015"), "num_samples"] <-length(unique(suppl_files[["Darmanis_2015"]]$`Experiment name`))
df_meta[which(df_meta$projects=="Eze_2021"), "cell_type"] <- paste(unique(suppl_files[["Eze_2021"]]$`Cell Type`), collapse = ", ")
df_meta[which(df_meta$projects=="Han_2020"), "age"] <- paste(unique(suppl_files[["Han_2020"]]$Age), collapse = ", ")
df_meta[which(df_meta$projects=="Han_2020"), "sex"] <- paste(unique(suppl_files[["Han_2020"]]$Gender), collapse = ", ")
df_meta[which(df_meta$projects=="Muller_2018"), "age"] <- paste(unique(suppl_files[["Muller_2018"]]$Age), collapse = ", ")
df_meta[which(df_meta$projects=="Muller_2018"), "age"] <- paste("YEARS", df_meta[which(df_meta$projects=="Muller_2018"), "age"], sep = " ")
df_meta[which(df_meta$projects=="Muller_2018"), "sex"] <- paste(unique(suppl_files[["Muller_2018"]]$Gender), collapse = ", ")

write.csv(df_meta, "/Users/aurazelco/Desktop/Lund_MSc/Thesis/data/UCSC/meta_summary.csv")