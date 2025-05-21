#!/bin/bash

input_folder="/home/cyprich/data/raw_data" #Enter path to input folder
output_folder="/home/cyprich/data/artifact_filtered" #Enter path to output folder
summary_file="/home/cyprich/data/summary_files/summary.csv"

# Create output folder and summary file directories if they don't exist
mkdir -p "$output_folder"
mkdir -p "$(dirname "$summary_file")"

# Summary Header
echo "Input File,Output File,Total Lines,Filtered Lines,Deleted Lines" > "$summary_file"

# Array to store unique valid values
declare -A unique_values

# Loop through each .vcf.gz file in the input folder
for gz_file in "$input_folder"/*.vcf.gz; do
    # Extract the base filename without extension
    base_name=$(basename "$gz_file" .vcf.gz)
    # Temporary decompressed file
    decompressed_file="$input_folder/$base_name.vcf"

    # Decompress the .vcf.gz file
    echo "Decompressing $gz_file..."
    gunzip -c "$gz_file" > "$decompressed_file"
    echo "$gz_file decompressed to $decompressed_file"

    # Output File name
    output_file="$output_folder/${base_name}_out1.vcf"

    # Count the total number of data rows
    total_lines=$(grep -v "^#" "$decompressed_file" | wc -l)

    # Filter and save to output file
    echo "Starting filtering for $decompressed_file..."
    awk 'BEGIN {OFS="\t"}
        /^#/ {print; next}  # We keep header lines starting with '#'
        $1 ~ /^chr[1-9]$|^chr1[0-9]$|^chr2[0-2]$|^chrX$|^chrY$|^chrM$/ {print}' "$decompressed_file" > "$output_file"

    # Extract unique valid values from column 1
    while read -r line; do
        if [[ ! $line =~ ^# ]]; then
            col1=$(echo "$line" | awk '{print $1}')
            unique_values["$col1"]=1
        fi
    done < "$output_file"

    # Count the number of filtered data rows
    filtered_lines=$(grep -v "^#" "$output_file" | wc -l)
    # Calculate the number of removed lines
    lines_deleted=$((total_lines - filtered_lines))

    # Append the summary for the current file as a new row in the CSV
    echo "$decompressed_file,$output_file,$total_lines,$filtered_lines,$lines_deleted" >> "$summary_file"
    echo "Processed $decompressed_file"

    # Remove temporary decompressed file
    rm "$decompressed_file"
    echo "Removed temporary file $decompressed_file"
done

# Print the length of the array (number of unique values)
echo "Total unique valid values from column 1: ${#unique_values[@]}"
# Print the unique valid values from column 1
echo "Unique valid values from column 1:"
for value in "${!unique_values[@]}"; do
    echo "$value"
done

echo "Summary report saved to $summary_file"
