# Exploring the sex and gender differences in the human brain at single cell level

## Table of contents
* [BACKGROUND](#background)
* [PROGRAMMING LANGUAGES](#programming-languages)
* [DATASETS](#datasets)
  * [DISCO](#disco)
    * [Dataset source](#dataset-source)
    * [RDS preparation](#rds-preparation)
  * [UCSC](#ucsc)
    * [Dataset source](#dataset-source)
    * [RDS preparation](#rds-preparation)
* [DIFFERENTIAL EXPRESSION ANALYSIS](#differential-expression-analysis)
* [INTEGRATION ANALYSIS](#integration-analysis) 
* [SUPPLEMENTARY FIGURES](#supplementary-figures)

----------------------------------------------------------------------------------------------------------

## BACKGROUND

*Abstract from the paper*

**Background:** Sex and gender have only lately been systematically regarded as a biological variable in both pre-clinical and clinical investigations. The impact of sex and gender on a wide range of biological and psychological variables remains largely unknown during brain development and ageing disorders.

**Methods:** To systematically evaluate sex and gender (SG) differences at different development stages and ageing disorders at a single cell level, we gathered publicly available single-nucleus RNA-sequencing studies through human life span from second trimester of gestation until geriatric age in healthy individuals and from Alzheimer's disease (AD) and Multiple Sclerosis (MS) patients. In summary, we collected single cell data for a total of 419885 single nuclei from 161 human brain samples (72 females and 89 males) to identify and characterise SG-biased genes.

**Results:** We identified SG-biased genes in both females and males in 11 major brain cell types across 7 developmental stages and two brain disorders. SG-biased genes were located mostly on the autosomes, with an enrichment for Y-linked genes in males in some cell types. SG-biased genes in many cell types were enriched for cell type markers. Accordingly, SG-biased genes showed little overlap across cell types and developmental stages. Interestingly, there was extensive functional overlap across SG-biased genes in developmental stages. Female-biased genes were enriched for brain-related functions and processes, and male-biased genes were enriched for metabolic pathways. Common female-biased genes across cell types and developmental stages contained many mitochondrial genes. Investigation of hormonal influence identified thymosin targets enriched in male-biased genes, and SG-biased genes for both males and females were highly enriched for androgen (not oestrogen) response elements across cell types.

**Conclusion:** We systematically characterised SG differences in brain development and brain-related disorders at a single cell level, leading to the identification of hormonal influences likely establishing these differences as well as enriched pathways and functional categories likely contributing to the SG differences in brain-related disorders. We have further made the entire analysis available as a web resource for the scientific community by developing a web application which can be found here ([app](), [source](https://github.com/aurazelco/HumanBrainSexSingleCell)).

*Why SG-biased differential expressed genes?*

The reason why we chose the "SG-biased" and not "sex-biased" differential expressed genes (DEGs) terminology is because, while we compared samples based on the biological sex (F/M), we cannot known the exact cause of the differences. Such differences in expression could be due to either or both sex and gender, and therefore we used the "SG" terminology. 

----------------------------------------------------------------------------------------------------------

## PROGRAMMING LANGUAGES

The analysis was performed mainly in R, with some bash scripts and command line commands. As a general rule, the R scripts could technically be run from the command line, however we do not recommend it. This is mainly due to the fact that some analysis need some user input, e.g. the number of dimensions to create UMAPs in Seurat. Therefore, we recommend to check and run sections of R inside an IDE of choice (e.g. RStudio) for optimal results. 
   
----------------------------------------------------------------------------------------------------------

## DATASETS

Below descriptions for the datasets included in this analysis can be found. 

### DISCO

#### Dataset source

The [DISCO](DISCO/) folder contains the scripts used on the DISCO dataset brain v1.0, found [here](https://www.immunesinglecell.org/atlasList). This dataset contains multiple publicly available single-cell and single-nucleus RNA-seq data, all from human samples. The advantage of using DISCO instead of the original projects is a unified cell annotation and metadata, which in turn can notably speed up the analysis. 

#### RDS preparation

This [folder](DISCO/RDS_preparation) contains the scripts to filter the DISCO dataset brain v1.0 according to a selection of criteria listed in details in the respective [README file](DISCO/README.md). Since the dataset can be downloaded as RDS file, there was no need to build a SeuratObject from raw data. 

Briefly, three projects were included in the analysis: [GSE157827](https://www.pnas.org/doi/10.1073/pnas.2008762117), [GSE174367](https://www.nature.com/articles/s41588-021-00894-z) and [PRJNA544731](https://www.nature.com/articles/s41586-019-1404-z). All the three projects had enough female and male samples, were mostly simial rin age and also contained patient samples, from Alzheimer's disease (GSE157827, GSE174367) and multiple sclerosis (PRJNA544731). 


### UCSC

#### Dataset source

The [UCSC](UCSC/) folder contains the scripts used on the datasets retrieved from the [UCSC Cell Browser](https://cells-test.gi.ucsc.edu). Among all datasets present, we selected the following dataset: Velmeshev et al. 2022 (bioRXiv) - Single-cell analysis of prenatal and postnatal human cortical development ([paper](https://www.biorxiv.org/content/10.1101/2022.10.24.513555v1.full.pdf), [UCSC dataset](https://cells-test.gi.ucsc.edu/?ds=cortical-dev+all))

This dataset contained not only fetal samples from the second and third trimester, but also data from the first years of life all the way into adulthood (individuals older than 20 years old). 

#### RDS preparation

This [folder](UCSC/RDS_preparation) contains the scripts to build the RDS from the expression matrix and metadata. More details can be found in the respective [README file](UCSC/README.md).

Because the original expression matrix was very memory-consuming, the matrix was split in multiple datasets according to the age range (as found in the metadata), and then each age was analyzed spearately. 

----------------------------------------------------------------------------------------------------------

## DIFFERENTIAL EXPRESSION ANALYSIS

The scripts to characterize the SG-biased DEGs can be found in the respective folders in [DISCO](DISCO/DEGs_individual_projects_adjust_pval) and [UCSC](UCSC/DEGs_adjust_pval). The analysis pipelines contain the same steps, but the scripts differ slightly because of different files structures. 

The workflows contained the following steps:
- filter the cell types with less than 10 SG-biased DEGs with significant adjusted p-value (Bonferroni correction)
- map the SG-biased DEGs against the genome to obtain information about the chromosome location
- investigate the presence of Xpar1 and Xpar2 genes (not included in the paper)
- calculate the X and Y chromosome enrichment in the SG-biased DEGs through hyper-geometric distribution
- calculate the percentage of androgen and estrogen response element binding sites in the SG-biased DEGs
- calculate the percentage of SG-biased DEGs which are conserved in other primate species

Each of these steps was run for each project/age in separate scripts. More details can be found in the [DISCO](DISCO/README.md) and [UCSC](UCSC/README.md) README files. 


----------------------------------------------------------------------------------------------------------

## INTEGRATION ANALYSIS

In this [folder](Integration/DEGs), there are the scripts used to compare the DEG results from both DISCO and UCSC, generated beforehand with the scripts found here ([DISCO](DISCO/DEGs_individual_projects_adjust_pval) and [UCSC](UCSC/DEGs_adjust_pval)). A brief description of the scripts can be found in the respective [README file](Integration/README.md). 


Additionally, in this [folder](Integration/Functional_analysis), the scripts used to run a functional analysis on the SG-biased DEGs can be found. More details [here](Integration/README.md). 


----------------------------------------------------------------------------------------------------------

## SUPPLEMENTARY FIGURES

This [folder](Suppl_files) contains the scripts to generate the plot for the thesis report **NOT** generated during the main analysis. A brief description of the scripts can be found in the respective [README file](Suppl_files/README.md). 
