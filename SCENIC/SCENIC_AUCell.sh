#!/bin/bash

# sets the flags for the arguments in the script
while getopts i:c:a: flag
do
    case "${flag}" in
        # where the original CSV files are
        i) input_directory=${OPTARG};;
        # where the CTX outputs are
        c) ctx=${OPTARG};;
        # where to save the AUCell outputs
        a) auc=${OPTARG};;
    esac
done

# for each CSV file in the input directory
for i in $input_directory*.csv; do
  # removes path so only filename is saved - used for output names
  id=${i%.csv}
  id=${id##*/}
  # prints the filename, so we can keep track which file is being analyzed
  echo $id
  # sets the output names for the GRN and CTX softwares as string variables
  ctx_out="$ctx$id".csv
  auc_out="$auc$id".csv
  # runs AUCell on the file and saves it in the AUCell output directory; -t is because the CSV file is transposed
  # respect to what SCENIC wants
  pyscenic aucell $i $ctx_out --output $auc_out  --num_workers 20 -t
done

