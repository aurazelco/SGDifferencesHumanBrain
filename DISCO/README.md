# DISCO Analysis

## Table of contents
* [Brief description of the RDS_preparation scripts](#brief-description-of-the-rds_preparation-scripts)
* [Brief description of the DEGs scripts](#brief-description-of-the-degs-scripts)
* [Brief description of the SCENIC scripts](#brief-description-of-the-scenic-scripts)


## Brief description of the RDS_preparation scripts

* [00_filtering_disco.R](RDS_preparation/00_filtering_disco.R) - script to filter the DISCO dataset
* [00_filtering_disco_final.R](RDS_preparation/00_filtering_disco_final.R) - script to filter the DISCO dataset -> polished version
* [01_general.R](RDS_preparation/01_general.R) - script to produce first informative plots about number of cells per celltype per sex, divided per project

## Brief description of the DEGs scripts

### DEGs_common

* [01B_plot_num_genes_func.R](DEGs/01B_plot_num_genes_func.R) - script to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs/01C_num_chr_func.R) - script to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [01E_CellMarker_func.R](DEGs/01E_CellMarker_func.R) - script to calculate and plot the percentage of DEGs which are also markers (marker list retrieved from [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/))
* [02A_Fisher_func.R](DEGs/02A_Fisher_func.R) - script to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs/02B_ARE_ERE_func.R) - script to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M
* [02C_Conservation_func.R](DEGs/02C_Conservation_func.R) - script to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts_DISCO_common.R](DEGs/all_scripts_DISCO_common.R) - script sourcing and running the previous scripts, 01A-02C, on the combined DEGs from all projects in a certain disease


### DEGs_individual_projects

* [01B_plot_num_genes_func.R](DEGs/01B_plot_num_genes_func.R) - script to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs/01C_num_chr_func.R) - script to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [01E_CellMarker_func.R](DEGs/01E_CellMarker_func.R) - script to calculate and plot the percentage of DEGs which are also markers (marker list retrieved from [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/))
* [02A_Fisher_func.R](DEGs/02A_Fisher_func.R) - script to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs/02B_ARE_ERE_func.R) - script to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M, but in the individual projects instead of on the intersected DEGs
* [02C_Conservation_func.R](DEGs/02C_Conservation_func.R) - script to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts_DISCO_proj.R](DEGs/all_scripts_DISCO_proj.R) - script sourcing and running the previous scripts, 01A-02C, on the DEGs from the individual projects in a certain disease

## Brief description of the SCENIC scripts

We used SCENIC to run a Gene Regulatory Network (GRN) analysis, followed their pipeline which can be found in the [original paper](https://doi.org/10.1038/s41596-020-0336-2) or on their [GitHub repo](https://github.com/aertslab/SCENICprotocol). 

The scripts below were used to execute the SCENIC pipeline and plot the results:
* [SCENIC_analysis_DISCO.sh](SCENIC/SCENIC_analysis_DISCO.sh) - bash script which reads the .csv files in the input folder, generate GRNBoos2  .tsv results in a second folder and finally combines the original .csv and the .tsv file to predict regulons using cisTarget - more information [here](#scenic_analysis)
* [SCENIC_merge_DISCO.sh](SCENIC/SCENIC_merge_DISCO.sh) - bash script which takes the individual .csv files for each project and sex combination, and merges them, in order to obtain one .csv file per project and sex - more information [here](#scenic_merge)
* [check_SCENIC_results.R](SCENIC/check_SCENIC_results_DISCO.R) -  R script to check the results from the SCENIC, more specifically to plot the overlap among runs, the overlap with [scGRNom results](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00908-9), and to plot the expression of the transcription factors (TFs) and targets (TGs) in the cell types. 

Additionally, in this folder the [scenic requirements](SCENIC/scenic_requirements.txt) file can be found, which was extracted from the conda environemt in which SCENIC was run. 

## SCENIC_analysis.sh

**Arguments of the script**

The script takes the following arguments:
* **-i** - input directory, where all the expression matrices (as csv files; rows = genes, column = cells) are found
* **-g** - directory where to save the outputs from the GRNBoost2 step
* **-c** - directory where to save the outputs from the cisTarget step
* **-a** - directory where to save the outputs from the AUCell step
* **-e** - directory where all extra files needed to run GRNBoost2 and cisTarget are found

**Example**

An example of how to run the script can be found below:

```shell
./SCENIC_analysis_DISCO.sh -i 0_input_dfs -g 1_GRN -c 2_CTX -a 3_AUCell -e ../../extra_files/
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please make sure that the csv input files are in this exact format (rows = genes, column = cells); if they are transposed, the bash script needs to be modified -> remove the *-t* argument in both *pyscenic grn* and *ctx* commands in the for loop



## SCENIC_merge.sh

**Arguments of the script**

The script takes the following arguments:
* **-i** - input directory, where all the expression matrices (as csv files; rows = genes, column = cells) are found
* **-n** - number of randomly sampled data we have for each project-sex combination
* **-o** - directory where to save the merged outputs
* **-s** - the sex to find the files to merge
* **-o** - the project to find the files to merge

**Example**

An example of how to run the script can be found below:

```shell
./SCENIC_merge_DISCO.sh -i MS/0_input_dfs/sampled_100_cells/ -n 3 -o MS/0_input_dfs/sampled_100_cells_all/ -s F -p PRJNA544731
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please make sure that the csv input files are all in the same format (e.g. rows = genes, column = cells); the script does not check the layout of the files before merging them
