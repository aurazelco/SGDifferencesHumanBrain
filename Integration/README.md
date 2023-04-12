# Integration analysis

## Table of contents
* [Brief description of the DEG integration scripts](#brief-description-of-the-deg-integration-scripts)
* [Brief description of the Functional_analysis integration scripts](#brief-description-of-the-functional_analysis-integration-scripts)


## Brief description of the DEGs scripts

* [Compare_DEGs.R](DEGs/DEGs.R), [Compare_DEGs_func.R](DEGs/DEGs_func.R) - script to compare the DEG results from both DISCO and UCSC; generates presence heatmaps across all analyzed ages and disease conditions (whether the gene is found in the group or not), and the counts of how many genes are shared in how may age/condition groups. 
* [Compare_ARE.R](DEGs/Compare_ARE.R), [Compare_DEGs_func.R](DEGs/Compare_ARE_func.R) - script to compare the ARE sites percentages across the conditions, separing the sexes, for each ct


### Second trimester integration

#### Single-cell RNA sequencing datasets

From the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu), we used the following dataset for the second trimester integration:
1. Nowakowski et al. 2017 - Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex ([paper](https://www.science.org/doi/epdf/10.1126/science.aap8809), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortex-dev))
2. Eze et al. 2021 - Heterogeneity of Human Neuroepithelial Cells and Early Radial Glia ([paper](https://www.nature.com/articles/s41593-020-00794-1), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=early-brain))

These two datasets all contained fetal samples, with both female and male samples. However, since the studies were likely not designed with a sex comparison in mind, we could find little number of age-matching samples between females and males within the same dataset. Therefore, we decided to integrate the two datasets and group the samples according to the gestational trimester, instead of by gestational week. This strategy also allowed for better comparison with the results from the Velmeshev analysis. 


## Brief description of the Functional_analysis scripts

* [Enrichment.R](Functional_analysis/Enrichment.R), [Enrichment_func.R](Functional_analysis/Enrichment_func.R) - script to run functional enrichment across cell types within each dataset, and across datasets within each cell type. 
* [Enrichment_M_2nd_trim.R](Functional_analysis/Enrichment_M_2nd_trim.R), [Enrichment_M_2nd_trim_func.R](Functional_analysis/Enrichment_M_2nd_trim_func.R) - script to analyze the functional enrichment (same as above) of the shared genes across cell types in the second trimester male-biased DEGs.
* [Comparison_wth_ref_datasets.R](Functional_analysis/Comparison_wth_ref_datasets.R), [Comparison_wth_ref_datasets_func.R](Functional_analysis/Comparison_wth_ref_datasets_func.R) - script to compare the SG-DEGs with reference datasets, namely: [McKenzie et al. 2018](https://www.nature.com/articles/s41598-018-27293-5), [Chlamydas et al. 2022](https://link.springer.com/article/10.1007/s00109-022-02227-x) and [SFARI Gene database](https://gene.sfari.org/database/human-gene/). 


The functiona enrichment included: 
1. gene ontology biological processes
2. Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways
3. disease enrichment (DisGeNET, DisGeNET CURATED, disease ontology)
4. GWAS Catalog 2019 for SNP enrichment
5. DSigDB for drug enrichment
6. TRANSFAC and JASPAR PWMs for transcription factor binding sites

