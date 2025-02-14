library(Seurat)
library(stringr)
library(tidyverse)

# General variables

min_cells <- 100

downsampling_size <- c(3, 5, 10, 50, 100, 250, 500)

n_downsampling <- 10

# Velmeshev 3rd trimester

trim_3rd <- readRDS("UCSC/Seurat_UCSC/Velmeshev_2022_3rd_trimester.rds")
trim_3rd@meta.data$sex_ct <- paste(trim_3rd@meta.data$sex, trim_3rd@meta.data$cluster_final, sep="_")
Idents(trim_3rd) <- "sex_ct"


final_groups <- read.csv(paste0("UCSC/DEGs_adjust_pval/Velmeshev_2022_3rd_trimester/outputs/final_filt_", min_cells, ".csv"))
sexes <- unique(final_groups$sex)
sexes <- c("F"="Female", "M"="Male")

"Excitatory neurons" %in% final_groups$ct

f_exc <- subset(trim_3rd, sex_ct=="Female_Excitatory neurons")
m_exc <- subset(trim_3rd, sex_ct=="Male_Excitatory neurons")

rm(trim_3rd, final_groups)

f_exc@meta.data$downsampling <- sample(x = 1:n_downsampling, size = ncol(f_exc), replace = TRUE)
Idents(f_exc) <- "downsampling"


m_exc@meta.data$downsampling <- sample(x = 1:n_downsampling, size = ncol(m_exc), replace = TRUE)
Idents(m_exc) <- "downsampling"


for (size in downsampling_size) {
  out <- paste("Extra_figures/DEGs_downsampling/Velmeshev_3rd_trim/degs", size, sep = "_")
  dir.create(out, recursive = T, showWarnings = F)
  f_exc_downsampled <- subset(f_exc, downsample = size)
  m_exc_downsampled <- subset(m_exc, downsample = size)
  
  downsampled <- merge(f_exc_downsampled, m_exc_downsampled)
  
  downsampled@meta.data$sex_ct_downsampling <- paste(downsampled@meta.data$sex_ct, downsampled@meta.data$downsampling, sep="_")
  
  Idents(downsampled) <- "sex_ct_downsampling"
  
  for (k in sort(unique(downsampled@meta.data$downsampling))) {
    id1 <- paste("Female_Excitatory neurons", k, sep="_")
    id2 <- paste("Male_Excitatory neurons", k, sep="_")
    deg1 <- FindMarkers(downsampled, 
                        ident.1 = id1, 
                        ident.2 = id2,
                        logfc.threshold = 0.25,
                        min.pct = 0.1,
                        only.pos = TRUE)
    write.csv(deg1, file.path(out, paste0("F_downsampled_", k, ".csv")))
    deg2 <- FindMarkers(downsampled,
                        ident.1 = id2, 
                        ident.2 = id1,
                        logfc.threshold = 0.25,
                        min.pct = 0.1,
                        only.pos = TRUE)
    write.csv(deg2, file.path(out, paste0("M_downsampled_", k, ".csv")))
  }
}


rm(f_exc, f_exc_downsampled, m_exc, m_exc_downsampled, deg1, deg2)

velm_all_degs <- lapply(downsampling_size, function(size) {
  out <- paste("Extra_figures/DEGs_downsampling/Velmeshev_3rd_trim/degs", size, sep = "_")
  degs <- lapply(list.files(out), function(x) {read.csv(file.path(out, x))})
  names(degs) <- stringr::str_remove_all(list.files(out, full.names = F, pattern = ".csv", recursive = F),
                                         "_downsampled|.csv")
  
  degs <- bind_rows(degs, .id = "run")
  
  degs <- degs %>%
    filter(p_val_adj <= 0.05, 
           avg_log2FC>= log2(1.2)) %>%
    rename("gene" = X) 
  
  return(degs)
  
})

names(velm_all_degs) <- paste0("n_", downsampling_size)


velm_degs <- bind_rows(velm_all_degs, .id="size")

# n_3, n_5, n_10 not included because no DEGs left after filtering

velm_all_comps <- lapply(velm_degs %>% pull(size) %>% unique(), function(n_size) {
  pairwise_comp <- velm_degs %>%
    filter(size == n_size) %>%
    mutate(presence = 1) %>%
    select(run, gene, presence) %>%
    pivot_wider(names_from = run, values_from = presence, values_fill = 0)
  
  f_comp <- lapply(pairwise_comp %>% select(contains("F_")) %>% colnames(), function(x1) {
    comp <- sapply(pairwise_comp %>% select(contains("F_")) %>% select(-all_of(x1)) %>% colnames(), function(x2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(x1), all_of(x2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  
  m_comp <- lapply(pairwise_comp %>% select(contains("M_")) %>% colnames(), function(y1) {
    comp <- sapply(pairwise_comp %>% select(contains("M_")) %>% select(-all_of(y1)) %>% colnames(), function(y2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(y1), all_of(y2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  f_comp <- bind_rows(f_comp) %>% mutate(sex = "F")
  m_comp <- bind_rows(m_comp) %>% mutate(sex = "M")
  
  
  comps <- rbind(f_comp, m_comp) %>%
    pivot_longer(!sex, names_to = "size", values_to = "overlap")
  
  return(comps)
})
names(velm_all_comps) <- velm_degs %>% pull(size) %>% unique()


# DISCO GSE157827

disco_filt <- readRDS("DISCOv1.0/brainV1.0_all_FM_filt.rds")
disco_filt <- subset(disco_filt, subset = project_id == "GSE157827")

final_groups <- read.csv(paste0("DISCOv1.0/DEGs_common/outputs/final_filt_", min_cells, ".csv"))

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

exc_neu <- names(unified_annotation[which(unified_annotation=="Excitatory neurons")])

final_groups <- final_groups %>%
  mutate(ct_unified = ct) %>%
  filter(ct_unified %in% exc_neu,
         proj == "GSE157827") %>%
  mutate(ct = "Excitatory neurons") %>%
  group_by(sex, disease) %>%
  summarise(tot_counts = sum(count))

unified_annotation <- as.data.frame(unified_annotation) %>%
                        rownames_to_column(var = "ct")

disco_rows <- rownames(disco_filt@meta.data)
disco_filt@meta.data <- left_join(disco_filt@meta.data, unified_annotation, by = "ct")
rownames(disco_filt@meta.data) <- disco_rows

disco_filt <- subset(disco_filt, subset = unified_annotation == "Excitatory neurons")

disco_filt@meta.data$sex_disease <- paste(disco_filt@meta.data$gender, disco_filt@meta.data$disease, sep = "_")

norm <- subset(disco_filt, subset = sex_disease %in% c("F_Normal", "M_Normal"))
ad <- subset(disco_filt, subset = sex_disease %in% c("F_Alzheimer's disease", "M_Alzheimer's disease"))

rm(disco_filt, final_groups)

norm@meta.data$sex_ct <- paste(norm@meta.data$gender, norm@meta.data$unified_annotation, sep = "_")
ad@meta.data$sex_ct <- paste(ad@meta.data$gender, ad@meta.data$unified_annotation, sep = "_")

Idents(norm) <- "sex_ct"
Idents(ad) <- "sex_ct"


for (i in seq(1, n_downsampling, by=1)) {
  for (size in downsampling_size) {
    out <- paste("Extra_figures/DEGs_downsampling/GSE157827_Healthy/degs", size, sep = "_")
    dir.create(out, recursive = T, showWarnings = F)
    
    downsampled <- subset(norm, downsample = size)
    
    id1 <- "F_Excitatory neurons"
    id2 <- "M_Excitatory neurons"
    deg1 <- FindMarkers(downsampled, 
                          ident.1 = id1, 
                          ident.2 = id2,
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          only.pos = TRUE)
    write.csv(deg1, file.path(out, paste0("F_downsampled_", i, ".csv")))
    deg2 <- FindMarkers(downsampled,
                          ident.1 = id2, 
                          ident.2 = id1,
                          logfc.threshold = 0.25,
                          min.pct = 0.1,
                          only.pos = TRUE)
    write.csv(deg2, file.path(out, paste0("M_downsampled_", i, ".csv")))
  }
}

for (i in seq(1, n_downsampling, by=1)) {
  for (size in downsampling_size) {
    out <- paste("Extra_figures/DEGs_downsampling/GSE157827_AD/degs", size, sep = "_")
    dir.create(out, recursive = T, showWarnings = F)
    
    downsampled <- subset(ad, downsample = size)
    
    id1 <- "F_Excitatory neurons"
    id2 <- "M_Excitatory neurons"
    deg1 <- FindMarkers(downsampled, 
                        ident.1 = id1, 
                        ident.2 = id2,
                        logfc.threshold = 0.25,
                        min.pct = 0.1,
                        only.pos = TRUE)
    write.csv(deg1, file.path(out, paste0("F_downsampled_", i, ".csv")))
    deg2 <- FindMarkers(downsampled,
                        ident.1 = id2, 
                        ident.2 = id1,
                        logfc.threshold = 0.25,
                        min.pct = 0.1,
                        only.pos = TRUE)
    write.csv(deg2, file.path(out, paste0("M_downsampled_", i, ".csv")))
  }
}



rm(norm, ad, downsampled, deg1, deg2, disco_rows, unified_annotation)

# Import DEGs

## Velmeshev

velm_all_degs <- lapply(downsampling_size, function(size) {
  out <- paste("Extra_figures/DEGs_downsampling/Velmeshev_3rd_trim/degs", size, sep = "_")
  degs <- lapply(list.files(out), function(x) {read.csv(file.path(out, x))})
  names(degs) <- stringr::str_remove_all(list.files(out, full.names = F, pattern = ".csv", recursive = F),
                                         "_downsampled|.csv")
  
  degs <- bind_rows(degs, .id = "run")
  
  degs <- degs %>%
    filter(p_val_adj <= 0.05, 
           avg_log2FC>= log2(1.2)) %>%
    rename("gene" = X) 
  
  return(degs)
  
})

names(velm_all_degs) <- paste0("n_", downsampling_size)


velm_degs <- bind_rows(velm_all_degs, .id="size")

# n_3, n_5, n_10 not included because no DEGs left after filtering

velm_all_comps <- lapply(velm_degs %>% pull(size) %>% unique(), function(n_size) {
  pairwise_comp <- velm_degs %>%
    filter(size == n_size) %>%
    mutate(presence = 1) %>%
    select(run, gene, presence) %>%
    pivot_wider(names_from = run, values_from = presence, values_fill = 0)
  
  f_comp <- lapply(pairwise_comp %>% select(contains("F_")) %>% colnames(), function(x1) {
    comp <- sapply(pairwise_comp %>% select(contains("F_")) %>% select(-all_of(x1)) %>% colnames(), function(x2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(x1), all_of(x2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  
  m_comp <- lapply(pairwise_comp %>% select(contains("M_")) %>% colnames(), function(y1) {
    comp <- sapply(pairwise_comp %>% select(contains("M_")) %>% select(-all_of(y1)) %>% colnames(), function(y2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(y1), all_of(y2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  f_comp <- bind_rows(f_comp) %>% mutate(sex = "F")
  m_comp <- bind_rows(m_comp) %>% mutate(sex = "M")
  
  
  comps <- rbind(f_comp, m_comp) %>%
    pivot_longer(!sex, names_to = "size", values_to = "overlap")
  
  return(comps)
})
names(velm_all_comps) <- velm_degs %>% pull(size) %>% unique()

## GSE157827 Healthy 

healthy_all_degs <- lapply(downsampling_size, function(size) {
  out <- paste("Extra_figures/DEGs_downsampling/GSE157827_Healthy/degs", size, sep = "_")
  degs <- lapply(list.files(out), function(x) {read.csv(file.path(out, x))})
  names(degs) <- stringr::str_remove_all(list.files(out, full.names = F, pattern = ".csv", recursive = F),
                                         "_downsampled|.csv")
  
  degs <- bind_rows(degs, .id = "run")
  
  degs <- degs %>%
    filter(p_val_adj <= 0.05, 
           avg_log2FC>= log2(1.2)) %>%
    rename("gene" = X) 
  
  return(degs)
  
})

names(healthy_all_degs) <- paste0("n_", downsampling_size)


healthy_degs <- bind_rows(healthy_all_degs, .id="size")

# n_3, n_5, n_10 not included because no DEGs left after filtering

healthy_all_comps <- lapply(healthy_degs %>% pull(size) %>% unique(), function(n_size) {
  pairwise_comp <- healthy_degs %>%
    filter(size == n_size) %>%
    mutate(presence = 1) %>%
    select(run, gene, presence) %>%
    pivot_wider(names_from = run, values_from = presence, values_fill = 0)
  
  f_comp <- lapply(pairwise_comp %>% select(contains("F_")) %>% colnames(), function(x1) {
    comp <- sapply(pairwise_comp %>% select(contains("F_")) %>% select(-all_of(x1)) %>% colnames(), function(x2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(x1), all_of(x2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  
  m_comp <- lapply(pairwise_comp %>% select(contains("M_")) %>% colnames(), function(y1) {
    comp <- sapply(pairwise_comp %>% select(contains("M_")) %>% select(-all_of(y1)) %>% colnames(), function(y2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(y1), all_of(y2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  f_comp <- bind_rows(f_comp) %>% mutate(sex = "F")
  m_comp <- bind_rows(m_comp) %>% mutate(sex = "M")
  
  
  comps <- rbind(f_comp, m_comp) %>%
    pivot_longer(!sex, names_to = "size", values_to = "overlap")
  
  return(comps)
})
names(healthy_all_comps) <- healthy_degs %>% pull(size) %>% unique()

## GSE157827 AD

ad_all_degs <- lapply(downsampling_size, function(size) {
  out <- paste("Extra_figures/DEGs_downsampling/GSE157827_AD/degs", size, sep = "_")
  degs <- lapply(list.files(out), function(x) {read.csv(file.path(out, x))})
  names(degs) <- stringr::str_remove_all(list.files(out, full.names = F, pattern = ".csv", recursive = F),
                                         "_downsampled|.csv")
  
  degs <- bind_rows(degs, .id = "run")
  
  degs <- degs %>%
    filter(p_val_adj <= 0.05, 
           avg_log2FC>= log2(1.2)) %>%
    rename("gene" = X) 
  
  return(degs)
  
})

names(ad_all_degs) <- paste0("n_", downsampling_size)


ad_degs <- bind_rows(ad_all_degs, .id="size")

# n_3, n_5, n_10 not included because no DEGs left after filtering

ad_all_comps <- lapply(ad_degs %>% pull(size) %>% unique(), function(n_size) {
  pairwise_comp <- ad_degs %>%
    filter(size == n_size) %>%
    mutate(presence = 1) %>%
    select(run, gene, presence) %>%
    pivot_wider(names_from = run, values_from = presence, values_fill = 0)
  
  f_comp <- lapply(pairwise_comp %>% select(contains("F_")) %>% colnames(), function(x1) {
    comp <- sapply(pairwise_comp %>% select(contains("F_")) %>% select(-all_of(x1)) %>% colnames(), function(x2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(x1), all_of(x2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  
  m_comp <- lapply(pairwise_comp %>% select(contains("M_")) %>% colnames(), function(y1) {
    comp <- sapply(pairwise_comp %>% select(contains("M_")) %>% select(-all_of(y1)) %>% colnames(), function(y2) {
      pairwise_comp %>%
        mutate(common = rowSums(pairwise_comp %>% select(all_of(y1), all_of(y2)))) %>%
        filter(common == 2) %>%
        nrow()
    })
    comp <- as.data.frame(comp)
    colnames(comp) <- n_size
    return(comp)
  })
  
  f_comp <- bind_rows(f_comp) %>% mutate(sex = "F")
  m_comp <- bind_rows(m_comp) %>% mutate(sex = "M")
  
  
  comps <- rbind(f_comp, m_comp) %>%
    pivot_longer(!sex, names_to = "size", values_to = "overlap")
  
  return(comps)
})
names(ad_all_comps) <- ad_degs %>% pull(size) %>% unique()

# Merge everything


all_comps <- bind_rows(list("Velmeshev_3rd_trim" = bind_rows(velm_all_comps),
                            "GSE157827_Healthy" = bind_rows(healthy_all_comps),
                            "GSE157827_AD" = bind_rows(ad_all_comps)), .id = "group") %>%
              mutate(size = factor(size, levels = c("n_50", "n_100","n_250","n_500")),
                     group = factor(group, levels = c("Velmeshev_3rd_trim", "GSE157827_Healthy","GSE157827_AD"))) 

ggplot(all_comps, aes(size, overlap, fill = sex )) +
  geom_boxplot() +
  labs(x= "Downsampling size", y = "Number of pairwise common genes", fill="Sex") +
  facet_wrap(~ group, nrow = 3) + 
  theme_classic()

ggsave("Extra_figures/DEGs_downsampling/pairwise_comparisons.pdf", dpi = 300, 
       width = 10, height = 15, units = "cm")


# Only Velmeshev plot

velm_overlap <- bind_rows(velm_all_comps) %>%
                mutate(size = factor(size,  c("n_50", "n_100","n_250","n_500")))

ggplot(velm_overlap, aes(size, overlap, fill = sex )) +
  geom_boxplot() +
  labs(x= "Downsampling size", y = "Number of pairwise common genes", fill="Sex") +
  theme_classic()

ggsave("Extra_figures/DEGs_downsampling/Velmeshev_3rd_trim_pairwise_comparisons.pdf", dpi = 300, 
       width = 10, height = 8, units = "cm")


##### Overlap in same downsampling size

bind_rows(list("Velmeshev_3rd_trim" = bind_rows(velm_all_degs),
               "GSE157827_Healthy" = bind_rows(healthy_all_degs),
               "GSE157827_AD" = bind_rows(ad_all_degs)), .id = "group")
 
overlap <- bind_rows(list("Velmeshev_3rd_trim" = bind_rows(velm_all_degs, .id = "size"),
                          "GSE157827_Healthy" = bind_rows(healthy_all_degs, .id = "size"),
                          "GSE157827_AD" = bind_rows(ad_all_degs, .id = "size")), .id = "group") %>%
            filter(size == "n_500") %>%
            separate(run, sep = "_", into = c("sex", "k_downsampling")) %>%
            group_by(sex, gene, group) %>%
            summarise(count = n()) %>%
            filter(count==10) %>%
            arrange(group, sex) %>%
            select(-count)

write.csv(overlap, file = "Extra_figures/DEGs_downsampling/common_genes.csv")            


overall_overlap <- overlap %>%
  group_by(sex, gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  filter(n > 1)
