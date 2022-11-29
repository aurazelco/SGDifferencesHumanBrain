#!/bin/bash

# sets the flags for the arguments in the script
while getopts i:n:s:o: flag
do
    case "${flag}" in
        # where the original CSV files are
        i) input_directory=${OPTARG};;
        # how many projects in folder
        n) num=${OPTARG};;
        # which sex to merge
        s) sex=${OPTARG};;
        # where to save the outputs
        o) out=${OPTARG};;
    esac
done

tmp=$(mktemp)
#files=("$input_directory$proj"_"$sex"_*_"$num_file".csv)

for ((num_file=1; num_file<=$num; num_file++)); do
  out_name="$out$sex"_"$num_file".csv
  pattern="$input_directory$sex"_*_"$num_file".csv
  files=( $pattern )
  id=${files[0]}
  cp $id $out_name
  for file in "${files[@]:1}"; do
    join -t ',' --nocheck-order $out_name  "$file" | column -t > "$tmp" && mv "$tmp" $out_name 
  done
  #paste -d "," "$input_directory$proj"_"$sex"_*_"$num_file".csv > $out_name
done

