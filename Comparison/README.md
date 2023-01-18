# Comparison analysis

## Table of contents
* [Brief description of the DEG comparison scripts](#brief-description-of-the-deg-comparison-scripts)
* [Brief description of the Functional_analysis comparison scripts](#brief-description-of-the-functional_analysis-comparison-scripts)


## Brief description of the DEGs scripts

* [Compare_DEGs.R](DEGs/Compare_DEGs.R), [Compare_DEGs_func.R](DEGs/Compare_DEGs_func.R) - script to compare the DEG results from both DISCO and UCSC; generates presence heatmaps across all analyzed ages and disease conditions (whether the gene is found in the group or not), and the counts of how many genes are shared in how may age/condition groups. 
* [Compare_ARE.R](DEGs/Compare_ARE.R), [Compare_DEGs_func.R](DEGs/Compare_ARE_func.R) - script to compare the ARE sites percentages across the conditions, separing the sexes, for each ct

## Brief description of the Functional_analysis scripts

* [Compare_Enrichment.R](Functional_analysis/Compare_Enrichment.R), [Compare_Enrichment_func.R](Functional_analysis/Compare_Enrichment_func.R) - script to analyze and compare the GO results: 1) for each ct and condition, F v M; 2) for each ct, F and M separately across conditions. 