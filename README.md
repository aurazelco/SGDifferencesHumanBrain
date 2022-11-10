# Repo of all main scripts generated during the MSc thesis

## Table of contents
* [DATASETS](#datasets)
  * [DISCO](#disco)
    * [DEGs analysis on DISCO dataset - F vs M](#degs-analysis-on-disco-dataset---f-vs-m)
* [SOFTWARES/ANALYSES](#softwares/analyses)
  * [Gene Regulatory Network Analysis](#gene-regulatory-network-analysis)
    * [SCENIC](#scenic)
* [THESIS DRAFT](#thesis-draft)

## LANGUAGES

The scripts are mainly in R, with some bash scripts and commands. As a general rule, the R script could technically be run from the command line, however we do not recommend it. This is mainly due to the fact that some analysis need some user input, e.g. the number of dimensions to create UMAPs in Seurat. Therefore, we recommend to check and run sections of R inside an IDE of choice (e.g. RStudio) for optimal results. 
   

## DATASETS

Below are listed the scripts for the analyses done for each individual dataset. 

### DISCO

The [DISCO](DISCO/) folder contains all scripts used on the DISCO dataset brain v1.0, found [here](https://www.immunesinglecell.org/atlasList). 

#### DEGs analysis on DISCO dataset - F vs M

This [folder](DISCO/DEGs) contains the analysis done on the DEGs between F and M started in August 2022. A brief description of the scripts can be found in the respective [README file](DISCO/DEGs/README.md). 

## SOFTWARES/ANALYSES

Below are listed general scripts done for all datasets. 

### Gene Regulatory Network Analysis

This analysis takes all genes expressed by the cells and search for regulatory networks

#### SCENIC

The [SCENIC](SCENIC/) folder contains a bash script to run the two main steps of SCENIC on all files contained in a folder, ran in October 2022. A brief description of the scripts can be found in the respective [README file](SCENIC/README.md). 

## THESIS DRAFT

This [folder](Thesis_draft) contains the scripts to generate the plot for the thesis report **NOT** generated during the analysis. A brief description of the scripts can be found in the respective [README file](Thesis_draft/README.md). 
