# Integration analysis

## Table of contents
* [Brief description of the DEG integration scripts](#brief-description-of-the-deg-integration-scripts)
* [Brief description of the Functional_analysis integration scripts](#brief-description-of-the-functional_analysis-integration-scripts)
* [Brief description of the Second_trimester integration scripts](#brief-description-of-the-second_trimester-integration-scripts)



## Brief description of the DEGs scripts

* [DEGs.R](DEGs/DEGs.R), [DEGs_func.R](DEGs/DEGs_func.R) - script to compare the DEG results from both DISCO and UCSC; generates presence heatmaps across all analyzed datasets (whether the gene is found in the group or not), and the counts of how many genes are shared in how many datasets. 
* [ARE_ERE.R](DEGs/ARE_ERE.R), [ARE_ERE_func.R](DEGs/ARE_ERE_func.R) - script to compare the ARE  and ERE sites percentages across the datasets, sex and cell types - integrates the results from the DEG analysis according to a unfiying annotation
* [Conservation.R](DEGs/Conservation.R), [Conservation_func.R](DEGs/Conservation_func.R) - script to compare the fractions of conserved SG-biased DEGs - integrates the results from the DEG analysis according to a unfiying annotation
* [DEGs_2nd_trimester.R](DEGs/DEGs_2nd_trimester.R), [DEGs_2nd_trimester_func.R](DEGs/DEGs_2nd_trimester_func.R) - script to compare the second trimster SG-biased DEGs from this analysis with SG-biased DEGs from bulk-RNA seq ([O'Brien et al. 2018](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1567-1)) and single-cell RNA seq (Eze et al. 2021, Nowakowski et al. 2017 - see below)
* [genes_location.R](DEGs/genes_location.R), [genes_location_func.R](DEGs/genes_location_func.R) - script to calculate the hyper-geometric enrichment of cellular compartments among the SG-biased DEGs using [Thul et al. 2017](https://www.science.org/doi/10.1126/science.aal3321)
* [hormones.R](DEGs/hormones.R), [hormones_func.R](DEGs/hormones_func.R) - script to calculate the hyper-geometric enrichment of hormonal targets among the SG-biased DEGs using Hormone-Gene version 1 from Jadhav et al. 2022 ([paper](https://academic.oup.com/bioinformatics/article/38/20/4771/6674503), [GitHub repo](https://github.com/BIRDSgroup/BioEmbedS))
* [shared_DEGs.R](DEGs/shared_DEGs.R), [shared_DEGs_func.R](DEGs/shared_DEGs_func.R) - script to generate presence heatmaps (whether the gene is present or not in the SG-biased DEGs) of the SG-biased DEGs shared by at least 75% of the cell type in each dataset, divided by chromosomal annotation (autosome, X and Y)
* [XY_Enrichment.R](DEGs/XY_Enrichment.R), [XY_Enrichment_func.R](DEGs/XY_Enrichment_func.R) - script to calculate the hyper-geometric enrichment of X and Y chromosome genes among the SG-biased DEGs - integrates the results from the DEG analysis according to a unfiying annotation


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


## Brief description of the Second_trimester scripts

* [DEGs_Eze_Nowa.R](Second_trimester/DEGs_Eze_Nowa.R) - script to generate the SG-biased DEGs used in the [DEGs_2nd_trimester.R](DEGs/DEGs_2nd_trimester.R) script