# Brief description of the scripts used to analyze datasets from UCSC Browser

* [00_metadata_parsing.R](00_metadata_parsing.R) - R script to collect metadata about the UCSC potential datasets to analyze
* [00_UCSC_Seurat.R](00_UCSC_Seurat.R) - R script to create the SeuratObjects from the expression matrices for the Nowakowski, Eze and van Bruggen projects
* [split_big_ExprMtx.sh](split_big_ExprMtx.sh) - bash script to split a heavy TSV.GZ file into multiple smaller TSV.GZ files according to indexes provided in .TXT files
* [00_UCSC_Seurat_Velmeshev.R](00_UCSC_Seurat_Velmeshev.R) - R script to create the SeuratObjects from the filtered expression matrix from Velmeshev

## Scripts

### split_big_ExprMtx

Before we ran this bash script for the Velmeshev dataset, we first removed all rows which contained only 0s, which represent genes with no expression across all cells. We did it by running the following command:

```shell
gunzip -c exprMatrix_Velmeshev_2022.tsv.gz | awk '{s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' | gzip > exprMtx_filt_Velmeshev_2022.tsv.gz
```

The following script does not mantain the header; however, we can easily retrieve it from the metadata in R, since the order of the samples in the column in the expression matrix and in the rows in the metadata *cell* column is the same. 
