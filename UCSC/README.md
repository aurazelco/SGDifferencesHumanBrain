# Brief description of the scripts used to analyze datasets from UCSC Browser

* [00_metadata_parsing.R](00_metadata_parsing.R) - R script to collect metadata about the UCSC potential datasets to analyze
* [00_UCSC_Seurat.R](00_UCSC_Seurat.R) - R script to create the SeuratObjects from the expression matrices for the Nowakowski, Eze and van Bruggen projects
* [split_big_ExprMtx.sh](split_big_ExprMtx.sh) - bash script to split a heavy TSV.GZ file into multiple smaller TSV.GZ files according to indexes provided in .TXT files
* [00_UCSC_Seurat_Velmeshev.R](00_UCSC_Seurat_Velmeshev.R) - R script to create the SeuratObjects from the filtered expression matrix from Velmeshev
* [00_UCSC_integration.R](00_UCSC_integration.R) -  R script to integrate Eze_2021 and Nowakowski_2017 datasets

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
