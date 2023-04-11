# DISCO Analysis

## Table of contents
* [Brief description of the RDS_preparation scripts](#brief-description-of-the-rds_preparation-scripts)
* [Brief description of the DEGs scripts](#brief-description-of-the-degs-scripts)


## Brief description of the RDS preparation scripts

* [00_filtering_disco.R](RDS_preparation/00_filtering_disco.R) - script to filter the DISCO dataset
* [00_filtering_disco_final.R](RDS_preparation/00_filtering_disco_final.R) - script to filter the DISCO dataset -> polished version
* [01_general.R](RDS_preparation/01_general.R) - script to produce first informative plots about number of cells per celltype per sex, divided per project

## Brief description of the DEGs scripts

* [01B_plot_num_genes_func.R](DEGs_individual_projects_adjust_pval/01B_plot_num_genes_func.R) - script to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs_individual_projects_adjust_pval/01C_num_chr_func.R) - script to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs_individual_projects_adjust_pval/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [01E_CellMarker_func.R](DEGs_individual_projects_adjust_pval/01E_CellMarker_func.R) - script to calculate and plot the percentage of DEGs which are also markers (marker list retrieved from [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/))
* [02A_Fisher_func.R](DEGs_individual_projects_adjust_pval/02A_Fisher_func.R) - script to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs_individual_projects_adjust_pval/02B_ARE_ERE_func.R) - script to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M, but in the individual projects instead of on the intersected DEGs
* [02C_Conservation_func.R](DEGs_individual_projects_adjust_pval/02C_Conservation_func.R) - script to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts_DISCO_proj.R](DEGs_individual_projects_adjust_pval/all_scripts_DISCO_proj.R) - script sourcing and running the previous scripts, 01A-02C, on the DEGs from the individual projects in a certain disease

