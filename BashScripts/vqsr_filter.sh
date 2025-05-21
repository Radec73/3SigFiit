#!/bin/bash

input_folder="/home/cyprich/data/VQSR_Processed"
output_folder="/home/cyprich/data/VQSR_Filtered"

for input_file in "$input_folder"/*.vcf; do

    file_name=$(basename "$input_file" .vcf)
    output_file="${output_folder}/${file_name}_vqsr_filtered.vcf"

    # Filter the VCF: keep only rows where 7th column is "PASS"
    awk 'BEGIN {OFS="\t"} !/^#/ {if ($7 == "PASS") print $0} /^#/ {print $0}' "$input_file" > "$output_file"

    echo "Filtered VCF saved to: $output_file"
done

echo "Processing complete.