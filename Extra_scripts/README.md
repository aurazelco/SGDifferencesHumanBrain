# DISCO Analysis

## Table of contents
* [Brief description of the Cistrome_DB scripts](#brief-description-of-the-cistrome_db-scripts)
* [Brief description of the Pybiomart scripts](#brief-description-of-the-pybiomart-scripts)



## Brief description of the Cistrome_DB scripts

* [cistrome_DB_TFs.R](Cistrome_DB/cistrome_DB_TFs.R), [cistrome_DB_TFs_func.R](istrome_DB/cistrome_DB_TFs_func.R)  - scripts to find the overlapping TFs between the GRNBoost2 SCENIC results (for both DISCO and UCSC) and the [Cistrome DB](http://cistrome.org/db/#/)

The file containing the TFs was copy-pasted from the Cistrome website (http://cistrome.org/db/#/) on 2023/01/11, selecting *Homo sapiens* as species and searching for the keyowrd *brain*. 

## Brief description of the Pybiomart scripts
* [02C_ENSEMBL_Biomart.ipynb](Pybiomart/02C_ENSEMBL_Biomart.ipynb) - jupyter notebook to create the mammalian reference dataset used in the DEG analysis. 

The jupyter notebook was ran inside a conda environment, of which the requirements can be found [here](../Conda_environments/Pybiomart/requirements.txt)


