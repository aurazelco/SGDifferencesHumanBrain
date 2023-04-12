# UCSC Analysis

## Table of contents
* [Brief description of the scripts used to prepare the RDS from UCSC](#brief-description-of-the-rds_preparation-scripts)
	* [Bash scripts](#bash-scripts)
		* [Bash command to filter the expression matrix](#bash-command-to-filter-the-expression-matrix)
		* [split_big_ExprMtx](#split_big_exprmtx)
	* [Dataset sources for second trimester integration](#dataset-sources-for-second-trimester-integration)
* [Brief description of the DEGs scripts](#brief-description-of-the-degs-scripts)
* [Brief description of the SCENIC scripts](#brief-description-of-the-scenic-scripts)


## Brief description of the scripts used to prepare the RDS from UCSC

* [Eze_Nowa_integration.R](RDS_preparation/Eze_Nowa_integration.R) -  R script to integrate Eze_2021 and Nowakowski_2017 datasets
* [Eze_Nowa_rds.R](RDS_preparationEze_Nowa_rds.R) - R script to create the SeuratObjects from the expression matrices for the Nowakowski and Eze projects
* [split_big_ExprMtx.sh](RDS_preparation/split_big_ExprMtx.sh) - bash script to split a heavy TSV.GZ file into multiple smaller TSV.GZ files according to indexes provided in .TXT files
* [UCSC_integration_2nd_trimester_all.R](RDS_preparation/UCSC_integration_2nd_trimester_all.R) -  R script to integrate all datasets with data from the 2nd trimester (Eze, Nowakowski, Velmeshev) - just to verify that we do not have separation based on dataset alone
* [UCSC_metadata_parsing.R](RDS_preparation/UCSC_metadata_parsing.R) - R script to collect metadata about the UCSC potential datasets to analyze
* [UCSC_Velmeshev_all_ages_furu.R](RDS_preparation/UCSC_Velmeshev_all_ages_furu.R) - R script to integrate all previous RDS SeuratObjects from Velmeshev
* [Velmeshev_furu_rds.R](RDS_preparation/Velmeshev_furu_rds.R) -  R script to create the RDS SeuratObjects for the larger Velmeshev subsets (2nd trimester, 10-20 years)
* [Velmeshev_rds.R](RDS_preparation/Velmeshev_rds.R) -  R script to create the RDS SeuratObjects for the smaller Velmeshev subsets, also prepares files for [split_big_ExprMtx.sh](RDS_preparation/split_big_ExprMtx.sh)

The RDS from UCSC were filtered according to the following steps:
- excluded projects that did not contain data from both sexes (at least 2-3/samples per sex) -> Velmeshev 4-10 years excluded, and Eze-Nowakowski integrated
- for each dataset, kept only the cell types that had samples from both sexes
- lastly, cell types with less than 100 cells per project and sex were not analyzed


### Bash scripts

#### Bash command to filter the expression matrix

Before we ran this bash script for the Velmeshev dataset, we first removed all rows which contained only 0s, which represent genes with no expression across all cells. We did it by running the following command:

```shell
gunzip -c exprMatrix_Velmeshev_2022.tsv.gz | awk '{s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' | gzip > exprMtx_filt_Velmeshev_2022.tsv.gz
```

The following script does not mantain the header; however, we can easily retrieve it from the metadata in R, since the order of the samples in the column in the expression matrix and in the rows in the metadata *cell* column is the same. 


#### split_big_ExprMtx

**Arguments of the script**

The script takes the following arguments:
* **-i** - input expression matrix to be split in multiple subsets (extract different columns)
* **-b** - input directory where the indexes of the columns to be extracted for each subset are saved (separate .TXT files, indexes listed in one line with commas as separator -> see  the object split_indexes in [00_UCSC_Seurat_Velmeshev.R](00_UCSC_Seurat_Velmeshev.R))
* **-o** - the output directory where to save the split expression matrices

**Example**

An example of how to run the script can be found below:

```shell
./split_big_ExprMtx.sh -i exprMtx_filt_Velmeshev_2022.tsv.gz -b Velmeshev_split -o Velmeshev_outs
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please note that the script may raise errors if the indexes are not sorted in ascending order
* In our case, we had to split the indexes of three subsets in two since the list of columns was too long 
* Please note that for 3 subsets (0_1_years_1, 2nd_trimester_2, 3rd_trimester_2), the expression matrices were obtained in R instead of bash, since the script was still raising errors - we ackowledge this is a limitation of the bash script and if this step becomes a routine step in the workflow, improvements in the script will be tested


### Dataset sources for second trimester integration

From the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu), we used the following dataset for the second trimester integration:
1. Nowakowski et al. 2017 - Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex ([paper](https://www.science.org/doi/epdf/10.1126/science.aap8809), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortex-dev))
2. Eze et al. 2021 - Heterogeneity of Human Neuroepithelial Cells and Early Radial Glia ([paper](https://www.nature.com/articles/s41593-020-00794-1), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=early-brain))

These two datasets all contained fetal samples, with both female and male samples. However, since the studies were likely not designed with a sex comparison in mind, we could find little number of age-matching samples between females and males within the same dataset. Therefore, we decided to integrate the two datasets and group the samples according to the gestational trimester, instead of by gestational week. This strategy also allowed for better comparison with the results from the Velmeshev analysis. The scripts for the RDS and integration are found above. 


## Brief description of the DEGs scripts

* [01B_plot_num_genes_func.R](DEGs_individual_projects_adjust_pval/01B_plot_num_genes_func.R) - functions to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs_individual_projects_adjust_pval/01C_num_chr_func.R) - functions to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs_individual_projects_adjust_pval/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [02A_Fisher_func.R](DEGs_individual_projects_adjust_pval/02A_Fisher_func.R) - functions to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs_individual_projects_adjust_pval/02B_ARE_ERE_func.R) - functions to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M, but in the individual projects instead of on the intersected DEGs
* [02C_Conservation_func.R](DEGs_individual_projects_adjust_pval/02C_Conservation_func.R) - functions to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts_Velmeshev_2nd_trim.R](DEGs/all_scripts_Velmeshev_2nd_trim.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 2nd trimester
* [all_scripts_Velmeshev_3rd_trim.R](DEGs/all_scripts_Velmeshev_3rd_trim.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 3rd trimester
* [all_scripts_Velmeshev_0_1_years.R](DEGs/all_scripts_Velmeshev_0_1_years.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 0-1 years
* [all_scripts_Velmeshev_1_2_years.R](DEGs/all_scripts_Velmeshev_1_2_years.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 1-2 years
* [all_scripts_Velmeshev_2_4_years.R](DEGs/all_scripts_Velmeshev_2_4_years.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 2-4 years
* [all_scripts_Velmeshev_10_20_years.R](DEGs/all_scripts_Velmeshev_10_20_years.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev 10-20 years
* [all_scripts_Velmeshev_Adult.R](DEGs/all_scripts_Velmeshev_Adult.R) - script containing the previous scripts, 01A-02C, to be run for the Velmeshev Adults
* [DEGs_Eze_Nowa.R](Second_trimester/DEGs_Eze_Nowa.R) - script to generate the SG-biased DEGs used in the [DEGs_2nd_trimester.R](DEGs/DEGs_2nd_trimester.R) script




