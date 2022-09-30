# Repo of all main scripts generated during the MSc thesis

## Table of contents
* [DISCO](#disco)
  * [DEGs analysis on DISCO dataset - F vs M](#degs-analysis-on-disco-dataset---f-vs-m)


## DISCO

The [DISCO](DISCO/) folder contains all scripts used on the DISCO dataset brain v1.0. 

### DEGs analysis on DISCO dataset - F vs M

This [folder](DISCO/DEGs) contains the analysis done on the DEGs between F and M started in August 2022.  

Brief description of the scripts:  
* [00_filtering_disco.R](DISCO/DEGs/00_filtering_disco.R) - script  to filter the DISCO dataset
* [01_general.R](DISCO/DEGs/01_general.R) - script  to produce first informative plots about number of cells per celltype per sex, divided per project
* [01A_generate_DEGs.R](DISCO/DEGs/01A_generate_DEGs.R) - script  to generate all the DEGs, one for F and one for M DEGs in each subtype and project
* [01B_plot_num_genes.R](DISCO/DEGs/01B_plot_num_genes.R), [01B_plot_num_genes_func.R](DISCO/DEGs/01B_plot_num_genes_func.R) - script to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr.R](DISCO/DEGs/01C_num_chr.R), [01C_num_chr_func.R](DISCO/DEGs/01C_num_chr_func.R) - script to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2.R](DISCO/DEGs/01D_Xpar1,2.R), [01D_Xpar1,2_func.R](DISCO/DEGs/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [02A_Fisher.R](DISCO/DEGs/02A_Fisher.R), [02A_Fisher_func.R](DISCO/DEGs/02A_Fisher_func.R) - script to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE.R](DISCO/DEGs/02B_ARE_ERE.R), [02B_ARE_ERE_func.R](DISCO/DEGs/02B_ARE_ERE_func.R) - script to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M
* [02C_Conservation.R](DISCO/DEGs/02C_Conservation.R), [02C_Conservation_func.R](DISCO/DEGs/02C_Conservation_func.R) - script to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts.R](DISCO/DEGs/all_scripts.R) - script containing the previous scripts, 01A-02C
