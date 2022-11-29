# Brief description of the scripts used to prepare the RDS from UCSC

* [00_metadata_parsing.R](00_metadata_parsing.R) - R script to collect metadata about the UCSC potential datasets to analyze
* [00_UCSC_Seurat.R](00_UCSC_Seurat.R) - R script to create the SeuratObjects from the expression matrices for the Nowakowski, Eze and van Bruggen projects
* [split_big_ExprMtx.sh](split_big_ExprMtx.sh) - bash script to split a heavy TSV.GZ file into multiple smaller TSV.GZ files according to indexes provided in .TXT files
* [00_UCSC_Seurat_Velmeshev.R](00_UCSC_Seurat_Velmeshev.R) - R script to create the SeuratObjects from the filtered expression matrix from Velmeshev
* [00_UCSC_integration.R](00_UCSC_integration.R) -  R script to integrate Eze_2021 and Nowakowski_2017 datasets
* [00_UCSC_Seurat_Velmeshev_furu.R](00_UCSC_Seurat_Velmeshev_furu.R) -  R script to create the RDS SeuratObjects for the larger Velmeshev subsets (2nd trimester, 10-20 years)
* [00_UCSC_Velmeshev_integration_furu.R](00_UCSC_Velmeshev_integration_furu.R) - R script to integrate all previous RDS SeuratObjects from Velmeshev
* [00_UCSC_integration_2nd_trimester_all.R](00_UCSC_integration_2nd_trimester_all.R) -  R script to integrate all datasets with data from the 2nd trimester (Eze, Nowakowski, Velmeshev)
* [00_UCSC_Velmeshev_analysis.R](00_UCSC_Velmeshev_analysis.R) -  R script to analyze the Velmeshev subsets        
* [00_UCSC_overall_plots.R](00_UCSC_overall_plots.R) -  R script to plot information about the Velmeshev subsets                


## Scripts

### Bash command to filter the expression matrix

Before we ran this bash script for the Velmeshev dataset, we first removed all rows which contained only 0s, which represent genes with no expression across all cells. We did it by running the following command:

```shell
gunzip -c exprMatrix_Velmeshev_2022.tsv.gz | awk '{s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' | gzip > exprMtx_filt_Velmeshev_2022.tsv.gz
```

The following script does not mantain the header; however, we can easily retrieve it from the metadata in R, since the order of the samples in the column in the expression matrix and in the rows in the metadata *cell* column is the same. 



### Bash - split_big_ExprMtx

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


# Brief description of the DEGs scripts

* [01A_generate_DEGs.R](DEGs/01A_generate_DEGs.R) - script  to generate all the DEGs, one for F and one for M DEGs in each subtype and project
* [01B_plot_num_genes_func.R](DEGs/01B_plot_num_genes_func.R) - script to extract the common DEGs among F and M DEGs per each subtype, and plot the results
* [01C_num_chr_func.R](DEGs/01C_num_chr_func.R) - script to map the DEGs obtained from 01B to the genome and plot the fraction of DEGs belonging to X, Y or autosomial chromosome; it also plots a heatmap of the X- and Y-genes expression across the different subtypes
* [01D_Xpar1,2_func.R](DEGs/01D_Xpar1,2_func.R) - script to calculate and plot the number of genes belonging to Xpar1 and Xpar2
* [01E_CellMarker_func.R](DEGs/01E_CellMarker_func.R) - script to calculate and plot the percentage of DEGs which are also markers (marker list retrieved from [CellMarker](http://bio-bigdata.hrbmu.edu.cn/CellMarker/))
* [02A_Fisher_func.R](DEGs/02A_Fisher_func.R) - script to calculate the enrichment of sex chromosomes compared to autosomial genes - used in 01C to add significance to plot
* [02B_ARE_ERE_func.R](DEGs/02B_ARE_ERE_func.R) - script to calculate and plot the ARE and ERE sites percentages in the cell types, separated F and M
* [02C_Conservation_func.R](DEGs/02C_Conservation_func.R) - script to plot the conserved fraction of DEGs across mammals and primates
* [all_scripts.R](DEGs/all_scripts.R) - script containing the previous scripts, 01A-02C

# Brief description of the SCENIC scripts

We used SCENIC to run a Gene Regulatory Network (GRN) analysis, followed their pipeline which can be found in the [original paper](https://doi.org/10.1038/s41596-020-0336-2) or on their [GitHub repo](https://github.com/aertslab/SCENICprotocol). 

The scripts below were used to execute the SCENIC pipeline and plot the results:
* [SCENIC_analysis_Eze_Nowa.sh](SCENIC/SCENIC_analysis_Eze_Nowa.sh) - bash script which reads the .csv files in the input folder, generate GRNBoos2  .tsv results in a second folder and finally combines the original .csv and the .tsv file to predict regulons using cisTarget - more information [here](#scenic_analysis)
* [SCENIC_merge_Eze_Nowa.sh](SCENIC/SCENIC_merge_Eze_Nowa.sh) - bash script which takes the individual .csv files for each project and sex combination, and merges them, in order to obtain one .csv file per project and sex - more information [here](#scenic_merge)
* [check_SCENIC_Eze_Nowa_results.R](SCENIC/check_SCENIC_Eze_Nowa_results.R), [check_SCENIC_results_func.R](SCENIC/check_SCENIC_results_func.R) -  R script to check the results from the SCENIC, more specifically to plot the overlap among runs, the overlap with [scGRNom results](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00908-9), and to plot the expression of the transcription factors (TFs) and targets (TGs) in the cell types. 


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
./SCENIC_analysis_Eze_Nowa.sh -i 0_input_dfs -g 1_GRN -c 2_CTX -a 3_AUCell -e ../../extra_files/
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
./SCENIC_merge_Eze_Nowa.sh -i MS/0_input_dfs/sampled_100_cells/ -n 3 -o MS/0_input_dfs/sampled_100_cells_all/ -s F -p PRJNA544731
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please make sure that the csv input files are all in the same format (e.g. rows = genes, column = cells); the script does not check the layout of the files before merging them


