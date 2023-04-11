# Repo of all main scripts generated during the MSc thesis

## Table of contents
* [BACKGROUND](#background)
* [PROGRAMMING LANGUAGES](#programming-languages)
* [DATASETS](#datasets)
  * [DISCO](#disco)
    * [Datasets source](#datasets-source)
    * [DEGs analysis on DISCO dataset](#degs-analysis-on-disco-dataset)
  * [UCSC](#ucsc)
    * [Datasets source](#datasets-source)
    * [DEGs analysis on UCSC datasets](#degs-analysis-on-ucsc-datasets)
* [INTEGRATION ANALYSIS](#integration-analysis) 
  * [DEGs](#degs)
    * [Second trimester integration](#second-trimester-integration)
    * []()
    * []()
    * []()
    * []()
    * []()
    * []()
    * []()
  * [Functional Enrichment](#functional-enrichment)
* [SUPPLEMENTARY FIGURES](#supplementary-figures)

----------------------------------------------------------------------------------------------------------

## BACKGROUND

*Abstract from the paper*

**Background:** Sex and gender have only lately been systematically regarded as a biological variable in both pre-clinical and clinical investigations. The impact of sex and gender on a wide range of biological and psychological variables remains largely unknown during brain development and ageing disorders.

**Methods:** To systematically evaluate sex and gender (SG) differences at different development stages and ageing disorders at a single cell level, we gathered publicly available single-nucleus RNA-sequencing studies through human life span from second trimester of gestation until geriatric age in healthy individuals and from Alzheimer's disease (AD) and Multiple Sclerosis (MS) patients. In summary, we collected single cell data for a total of 419885 single nuclei from 161 human brain samples (72 females and 89 males) to identify and characterise SG-biased genes.

**Results:** We identified SG-biased genes in both females and males in 11 major brain cell types across 7 developmental stages and two brain disorders. SG-biased genes were located mostly on the autosomes, with an enrichment for Y-linked genes in males in some cell types. SG-biased genes in many cell types were enriched for cell type markers. Accordingly, SG-biased genes showed little overlap across cell types and developmental stages. Interestingly, there was extensive functional overlap across SG-biased genes in developmental stages. Female-biased genes were enriched for brain-related functions and processes, and male-biased genes were enriched for metabolic pathways. Common female-biased genes across cell types and developmental stages contained many mitochondrial genes. Investigation of hormonal influence identified thymosin targets enriched in male-biased genes, and SG-biased genes for both males and females were highly enriched for androgen (not oestrogen) response elements across cell types.

**Conclusion:** We systematically characterised SG differences in brain development and brain-related disorders at a single cell level, leading to the identification of hormonal influences likely establishing these differences as well as enriched pathways and functional categories likely contributing to the SG differences in brain-related disorders. We have further made the entire analysis available as a web resource for the scientific community by developing a web application which can be found [here](https://www.immunesinglecell.org/atlasList).


----------------------------------------------------------------------------------------------------------

## LANGUAGES

The analysis was performed mainly in R, with some bash scripts and command line commands. As a general rule, the R scripts could technically be run from the command line, however we do not recommend it. This is mainly due to the fact that some analysis need some user input, e.g. the number of dimensions to create UMAPs in Seurat. Therefore, we recommend to check and run sections of R inside an IDE of choice (e.g. RStudio) for optimal results. 
   
----------------------------------------------------------------------------------------------------------

## DATASETS

Below descriptions for the datasets included in this analysis can be found. 

### DISCO

#### Dataset source

The [DISCO](DISCO/) folder contains the scripts used on the DISCO dataset brain v1.0, found [here](https://www.immunesinglecell.org/atlasList). 

#### DEGs analysis on DISCO dataset

This [folder](DISCO/DEGs) contains the analysis done on the DEGs between F and M. A brief description of the scripts can be found in the respective [README file](DISCO/README.md). 

----------------------------------------------------------------------------------------------------------

### UCSC

#### Dataset source

The [UCSC](UCSC/) folder contains the scripts used on the datasets retrieved from the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu). Among all datasets present, we selected the following dataset: Velmeshev et al. 2022 (bioRXiv) - Single-cell analysis of prenatal and postnatal human cortical development ([paper](https://www.biorxiv.org/content/10.1101/2022.10.24.513555v1.full.pdf), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortical-dev+all))

This dataset contained not only fetal samples from the second and third trimester, but also data from the first years of life all the way into adulthood (individuals older than 20 years old). 

#### DEGs analysis on UCSC datasets

This [folder](UCSC/DEGs) contains the analysis done on the DEGs between F and M. A brief description of the scripts can be found in the respective [README file](UCSC/README.md). 

----------------------------------------------------------------------------------------------------------

## INTEGRATION ANALYSIS

In this [folder](Integration/DEGs), there are the scripts used to compare the DEG results from both DISCO and UCSC, generated beforehand with the scripts found here ([DISCO](DISCO/DEGs) and [UCSC](UCSC/DEGs)). A brief description of the scripts can be found in the respective [README file](Integration/README.md). 


### Second trimester integration

#### Single-cell RNA sequencing datasets

From the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu), we used the following dataset for the second trimester integration:
1. Nowakowski et al. 2017 - Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex ([paper](https://www.science.org/doi/epdf/10.1126/science.aap8809), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortex-dev))
2. Eze et al. 2021 - Heterogeneity of Human Neuroepithelial Cells and Early Radial Glia ([paper](https://www.nature.com/articles/s41593-020-00794-1), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=early-brain))

These two datasets all contained fetal samples, with both female and male samples. However, since the studies were likely not designed with a sex comparison in mind, we could find little number of age-matching samples between females and males within the same dataset. Therefore, we decided to integrate the two datasets and group the samples according to the gestational trimester, instead of by gestational week. This strategy also allowed for better comparison with the results from the Velmeshev analysis. 


----------------------------------------------------------------------------------------------------------

## SUPPLEMENTARY FIGURES

This [folder](Suppl_files) contains the scripts to generate the plot for the thesis report **NOT** generated during the main analysis. A brief description of the scripts can be found in the respective [README file](Suppl_files/README.md). 
