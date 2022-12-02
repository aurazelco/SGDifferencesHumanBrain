# Repo of all main scripts generated during the MSc thesis

## Table of contents
* [BACKGROUND](#background)
* [PROGRAMMING LANGUAGES](#programming-languages)
* [DATASETS](#datasets)
  * [DISCO](#disco)
    * [Datasets source](#datasets-source)
    * [DEGs analysis on DISCO dataset](#degs-analysis-on-disco-dataset)
    * [Gene Regulatory Network Analysis - SCENIC](#gene-regulatory-network-analysis---scenic)
  * [UCSC](#ucsc)
    * [Datasets source](#datasets-source)
      * [Eze - Nowakowski](#eze---nowakowski)
      * [Velmeshev](#velmeshev)
    * [DEGs analysis on UCSC datasets](#degs-analysis-on-ucsc-datasets)
    * [Gene Regulatory Network Analysis - SCENIC](#gene-regulatory-network-analysis---scenic)
* [THESIS DRAFT](#thesis-draft)

----------------------------------------------------------------------------------------------------------

## BACKGROUND

## LANGUAGES

The scripts are mainly in R, with some bash scripts and commands. As a general rule, the R script could technically be run from the command line, however we do not recommend it. This is mainly due to the fact that some analysis need some user input, e.g. the number of dimensions to create UMAPs in Seurat. Therefore, we recommend to check and run sections of R inside an IDE of choice (e.g. RStudio) for optimal results. 
   

----------------------------------------------------------------------------------------------------------

## DATASETS

Below are listed the scripts for the analyses done for each individual dataset. 

### DISCO

#### Datasets source

The [DISCO](DISCO/) folder contains all scripts used on the DISCO dataset brain v1.0, found [here](https://www.immunesinglecell.org/atlasList). 

#### DEGs analysis on DISCO dataset

This [folder](DISCO/DEGs) contains the analysis done on the DEGs between F and M started in August 2022. A brief description of the scripts can be found in the respective [README file](DISCO/README.md). 

### Gene Regulatory Network Analysis

This analysis takes all genes expressed by the cells and search for regulatory networks

#### Gene Regulatory Network Analysis - SCENIC

The [SCENIC](DISCO/SCENIC/) folder contains a bash script to run the two main steps of SCENIC on all files contained in a folder, ran in October 2022. A brief description of the scripts can be found in the respective [README file](DISCO/README.md). 

----------------------------------------------------------------------------------------------------------

### UCSC

#### Datasets source

The [UCSC](UCSC/) folder contains all scripts used on the datasets retrieved from the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu). Among all datasets presetn, we selected the following three:
1. Nowakowski et al. 2017 - Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex ([paper](https://www.science.org/doi/epdf/10.1126/science.aap8809), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortex-dev))
2. Eze et al. 2021 - Heterogeneity of Human Neuroepithelial Cells and Early Radial Glia ([paper](https://www.nature.com/articles/s41593-020-00794-1), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=early-brain))
3. Velmeshev et al. 2022 (bioRXiv) - Single-cell analysis of prenatal and postnatal human cortical development ([paper](https://www.biorxiv.org/content/10.1101/2022.10.24.513555v1.full.pdf), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortical-dev+all))

###### Eze - Nowakowski

These datasets all contained fetal samples, with both female and male samples. However, since the studies were likely not designed with a sex comparison in mind, we could find little number of age-matching samples between females and males within the same dataset. Therefore, we decided to integrate the two datasets and group the samples according to the gestational trimester, instead of by gestational week. This strategy also allowed for better comparison with the results from the Velmeshev analysis. 

###### Velmeshev

This dataset contained not only fetal samples from the second and third trimester, but also data from the first years of life all the way into adulthood (individuals older than 20 years old). 

#### DEGs analysis on UCSC datasets

This [folder](UCSC/DEGs) contains the analysis done on the DEGs between F and M done in November/Decmber 2022. A brief description of the scripts can be found in the respective [README file](UCSC/README.md). 

#### Gene Regulatory Network Analysis - SCENIC

The [SCENIC](UCSC/SCENIC/) folder contains a bash script to run the two main steps of SCENIC on all files contained in a folder, ran in October 2022. A brief description of the scripts can be found in the respective [README file](UCSC/README.md). 

----------------------------------------------------------------------------------------------------------

## THESIS DRAFT

This [folder](Thesis_draft) contains the scripts to generate the plot for the thesis report **NOT** generated during the analysis. A brief description of the scripts can be found in the respective [README file](Thesis_draft/README.md). 
