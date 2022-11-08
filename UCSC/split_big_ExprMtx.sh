#!/bin/bash

# sets the flags for the arguments in the script
while getopts i:b:o: flag
do
    case "${flag}" in
        # where the original TSV.GZ file is
        i) input_TSV_file=${OPTARG};;
        # where the barcode indexes are
        b) bi_directory=${OPTARG};;
        # where to save the split TSV.GZ outputs
        o) out=${OPTARG};;
    esac
done

# temporary file, needed to generate the new ones as merged datasets
tmp=$(mktemp)


# for each CSV file in the input directory

for file in $bi_directory"/"*.txt; do
  # removes path so only filename is saved - used for output names
  id=${file%.txt}
  id=${id##*/}
  # prints the filename, so we can keep track which file is being analyzed
  echo $id
  # sets the output names for the GRN and CTX softwares as string variables
  #indexes=$(cat $file)
  indexes=$(<$file)
  #declare -i $indexes
  out_name="$out"/"$id".tsv
  echo $out_name
  gunzip -c $input_TSV_file | cut -f$indexes | gzip > $out_name
done

