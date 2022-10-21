#!/bin/bash

# sets the flags for the arguments in the script
while getopts i:n:s:p:o: flag
do
    case "${flag}" in
        # where the original CSV files are
        i) input_directory=${OPTARG};;
        # how many sampling we have for each project/sex combination
        n) num=${OPTARG};;
        # which sex to merge
        s) sex=${OPTARG};;
        # which project to merge
        p) proj=${OPTARG};;
        # where to save the outputs
        o) out=${OPTARG};;
    esac
done

# temporary file, needed to generate the new ones as merged datasets
tmp=$(mktemp)

# loops through the different versions of the datasets, e.g. num=3 means we have _1, _2 and _3 -> 3 different random sampling
for ((num_file=1; num_file<=$num; num_file++)); do
  # creates the output name
  out_name="$out$proj"_"$sex"_"$num_file".csv
  # extract the filenames matching the project, sex and sampling 
  pattern="$input_directory$proj"_"$sex"_*_"$num_file".csv
  files=( $pattern )
  # extract the first filename, since we still need the first file to contain the first column "Genes"
  id=${files[0]}
  # copy-stats the first file
  cp $id $out_name
  # for all files in the folder matching the pattern, excluded the first one that we have already copied
  for file in "${files[@]:1}"; do
    # merges the csv files without changing the order of the columns in the temporary file, and then copies this to the output name
    join -t ',' --nocheck-order $out_name  "$file" | column -t > "$tmp" && mv "$tmp" $out_name 
  done
done

