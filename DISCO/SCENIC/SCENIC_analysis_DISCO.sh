#!/bin/bash

# sets the flags for the arguments in the script
while getopts i:g:c:a:e: flag
do
    case "${flag}" in
        # where the original CSV files are
        i) input_directory=${OPTARG};;
        # where to save the GRN outputs
        g) grn=${OPTARG};;
        # where to save the CTX outputs
        c) ctx=${OPTARG};;
        # where to save the AUCell outputs
        a) auc=${OPTARG};;
        # where the extra files are located
        e) extra=${OPTARG};;
    esac
done

# defines as variables the path to the extra files - may raise errors if filename is changed
extra_hs="$extra"hs_hgnc_tfs.txt
extra_motif="$extra"motifs-v9-nr.hgnc-m0.001-o0.0.tbl
extra_genome="$extra"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

# for each CSV file in the input directory
for i in $input_directory*.csv; do
  # removes path so only filename is saved - used for output names
  id=${i%.csv}
  id=${id##*/}
  # prints the filename, so we can keep track which file is being analyzed
  echo $id
  # sets the output names for the GRN and CTX softwares as string variables
  grn_out="$grn$id".tsv
  ctx_out="$ctx$id".csv
  # runs GRNBoost2 on the file and saves it in the GRN output directory; -t is because the CSV file is transposed
  # respect to what SCENIC wants
  pyscenic grn --num_workers 20 --output $grn_out --method grnboost2 $i $extra_hs -t
  # runs CTX on the output from GRNBoost2; -t is because the raw CSV file is transposed
  pyscenic ctx $grn_out $extra_genome --annotations_fname $extra_motif --expression_mtx_fname $i --mode "dask_multiprocessing" --output $ctx_out --num_workers 20 --mask_dropouts -t
  auc_out="$auc$id".csv
  # runs AUCell on the file and saves it in the AUCell output directory; -t is because the CSV file is transposed
  # respect to what SCENIC wants
  pyscenic aucell $i $ctx_out --output $auc_out  --num_workers 20 -t
done

