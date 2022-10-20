# Brief description of the SCENIC script

We used SCENIC to run a Gene Regulatory Network (GRN) analysis, followed their pipeline which can be found in the [original paper](https://doi.org/10.1038/s41596-020-0336-2) or on their [GitHub repo](https://github.com/aertslab/SCENICprotocol). 

Due to the high number of files, we wrote a bash script to automate the SCENIC pipeline for hundreds of files. 

* [SCENIC_analysis.sh](SCENIC_analysis.sh) - bash script which reads the .csv files in the input folder, generate GRNBoos2  .tsv results in a second folder and finally combines the original .csv and the .tsv file to predict regulons using cisTarget. 

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
