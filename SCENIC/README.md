# Brief description of the SCENIC scripts

We used SCENIC to run a Gene Regulatory Network (GRN) analysis, followed their pipeline which can be found in the [original paper](https://doi.org/10.1038/s41596-020-0336-2) or on their [GitHub repo](https://github.com/aertslab/SCENICprotocol). 

The scripts below were used to execute the SCENIC pipeline and plot the results:
* [SCENIC_analysis.sh](SCENIC_analysis.sh) - bash script which reads the .csv files in the input folder, generate GRNBoos2  .tsv results in a second folder and finally combines the original .csv and the .tsv file to predict regulons using cisTarget - more information [here](#scenic_analysis)
* [SCENIC_merge.sh](SCENIC_merge.sh) - bash script which takes the individual .csv files for each project and sex combination, and merges them, in order to obtain one .csv file per project and sex - more information [here](#scenic_merge)
* [check_SCENIC_results.R](check_SCENIC_results.R), [check_SCENIC_results_func.R](check_SCENIC_results_func.R) -  R script to check the results from the SCENIC, more specifically to plot the overlap among runs, the overlap with [scGRNom results](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00908-9), and to plot the expression of the transcription factors (TFs) and targets (TGs) in the cell types. 


## SCENIC_analysis

**Arguments of the script**

The script takes the following arguments:
* **-i** - input directory, where all the expression matrices (as csv files; rows = genes, column = cells) are found
* **-g** - directory where to save the outputs from the GRNBoost2 step
* **-c** - directory where to save the outputs from the cisTarget step
* **-e** - directory where all extra files needed to run GRNBoost2 and cisTarget are found

**Example**

An example of how to run the script can be found below:

```shell
./SCENIC_analysis.sh -i 0_input_dfs -g 1_GRN -c 2_CTX -e ../../extra_files/
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please make sure that the csv input files are in this exact format (rows = genes, column = cells); if they are transposed, the bash script needs to be modified -> remove the *-t* argument in both *pyscenic grn* and *ctx* commands in the for loop



## SCENIC_merge

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
./SCENIC_merge.sh -i MS/0_input_dfs/sampled_100_cells/ -n 3 -o MS/0_input_dfs/sampled_100_cells_all/ -s F -p PRJNA544731
```

**Notes**

* Please note that if no argument is given, the script may raise error or overwrite files
* Please make sure that the csv input files are all in the same format (e.g. rows = genes, column = cells); the script does not check the layout of the files before merging them
