# Comparison analysis

## Table of contents
* [Brief description of the DEG comparison scripts](#brief-description-of-the-deg-comparison-scripts)
* [Brief description of the Functional_analysis comparison scripts](#brief-description-of-the-functional_analysis-comparison-scripts)


## Brief description of the DEGs scripts

* [Compare_DEGs.R](DEGs/Compare_DEGs.R), [Compare_DEGs_func.R](DEGs/Compare_DEGs_func.R) - script to compare the DEG results from both DISCO and UCSC; generates presence heatmaps across all analyzed ages and disease conditions (whether the gene is found in the group or not), and the counts of how many genes are shared in how may age/condition groups. 
* [Compare_ARE.R](DEGs/Compare_ARE.R), [Compare_DEGs_func.R](DEGs/Compare_ARE_func.R) - script to compare the ARE sites percentages across the conditions, separing the sexes, for each ct

## Brief description of the Functional_analysis scripts

* [Compare_Enrichment.R](Functional_analysis/Compare_Enrichment.R), [Compare_Enrichment_func.R](Functional_analysis/Compare_Enrichment_func.R) - script to analyze and compare the GO results: 1) for each ct and condition, F v M; 2) for each ct, F and M separately across conditions. 


### Second trimester integration

#### Single-cell RNA sequencing datasets

From the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu), we used the following dataset for the second trimester integration:
1. Nowakowski et al. 2017 - Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex ([paper](https://www.science.org/doi/epdf/10.1126/science.aap8809), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortex-dev))
2. Eze et al. 2021 - Heterogeneity of Human Neuroepithelial Cells and Early Radial Glia ([paper](https://www.nature.com/articles/s41593-020-00794-1), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=early-brain))

These two datasets all contained fetal samples, with both female and male samples. However, since the studies were likely not designed with a sex comparison in mind, we could find little number of age-matching samples between females and males within the same dataset. Therefore, we decided to integrate the two datasets and group the samples according to the gestational trimester, instead of by gestational week. This strategy also allowed for better comparison with the results from the Velmeshev analysis. 

