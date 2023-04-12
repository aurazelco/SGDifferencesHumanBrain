# DISCO Analysis

## Table of contents
* [Brief description of the RDS_preparation scripts](#brief-description-of-the-rds_preparation-scripts)
* [Brief description of the DEGs scripts](#brief-description-of-the-degs-scripts)
* [Workflow](#workflow)


## Brief description of the RDS preparation scripts

* [filtering_disco.R](RDS_preparation/filtering_disco.R) - script to filter the DISCO dataset

The original RDS from DISCO was filtered according to the following steps:
- excluded samples/projects that did not have sex metadata
- excluded projects that did not contain data from both sexes
- for each disease condition (normal/healthy, Alzheimer's disease, and multiple sclerosis) kept only the cell types that had samples from both sexes
- then excluded the cell types which had samples from only one project
- lastly, cell types with less than 100 cells per project and sex were not analyzed


## Brief description of the DEGs scripts

* [prep_DEG_files.R](DEGs_individual_projects_adjust_pval/prep_DEG_files.R) - script to prepare files for the DEG analysis, to be run **BEFORE** the DEGs analysis
* [01B_plot_num_genes_func.R](DEGs_individual_projects_adjust_pval/01B_plot_num_genes_func.R) - functions to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs_individual_projects_adjust_pval/01C_num_chr_func.R) - functions to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs_individual_projects_adjust_pval/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [01E_CellMarker_func.R](DEGs_individual_projects_adjust_pval/01E_CellMarker_func.R) - functions to calculate and plot the percentage of DEGs which are also markers (marker list retrieved from [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/))
* [02A_Fisher_func.R](DEGs_individual_projects_adjust_pval/02A_Fisher_func.R) - functions to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs_individual_projects_adjust_pval/02B_ARE_ERE_func.R) - functions to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M, but in the individual projects instead of on the intersected DEGs
* [02C_Conservation_func.R](DEGs_individual_projects_adjust_pval/02C_Conservation_func.R) - functions to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts_DISCO_proj_GSE157827.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_GSE157827.R) - script sourcing and running the previous scripts, 01A-02C, on the DEGs in an individual project
* [all_scripts_DISCO_proj_GSE174367.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_GSE174367.R) - script sourcing and running the previous scripts, 01A-02C, on the DEGs in an individual project
* [all_scripts_DISCO_proj_PRJNA544731.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_PRJNA544731.R) - script sourcing and running the previous scripts, 01A-02C, on the DEGs in an individual project


## Workflow

The scripts should be run in the following order:
1. [filtering_disco.R](RDS_preparation/filtering_disco.R)
2. [prep_DEG_files.R](DEGs_individual_projects_adjust_pval/prep_DEG_files.R)
3. [all_scripts_DISCO_proj_GSE157827.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_GSE157827.R), [all_scripts_DISCO_proj_GSE174367.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_GSE174367.R), [all_scripts_DISCO_proj_PRJNA544731.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj_PRJNA544731.R) (among these three, no order has to be followed)