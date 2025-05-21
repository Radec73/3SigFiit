#!/bin/bash

input_folder="/home/cyprich/data/af_filtered"
output_folder="/home/cyprich/data/final_data"
summary_file="/home/cyprich/data/final_data/summary.csv"

# Initialize the summary CSV file
echo "FileName,TotalLines,FilteredLines,LinesDeleted" > "$summary_file"
#REGEX naming convention
for tumor_sample in "$input_folder"/*_T[0-9]*-R[0-9]*_*.vcf; do

    #REGEX naming convention of one by one tumor sample
    patient_ID=$(echo "$tumor_sample" | sed -E 's/(.*_T_R[0-9]+).*/\1/')

    #Regex naming convention of its blood sample
    blood_sample="${patient_ID}_BS.vcf"

    filename=$(basename "$tumor_sample" .vcf)
    output_file="${output_folder}/${tumor_sample}_final.vcf"

    # Create a temporary file to store blood sample positions
    BLOOD_POSITIONS=$(mktemp)

    # Extract CHROM and POS from the blood VCF (skip headers)
    grep -v '^#' "$blood_sample" | awk '{print $1"\t"$2}' > "$BLOOD_POSITIONS"

    # Get the total number of data lines (non-header) in the input file
    total_lines=$(grep -v "^#" "$tumor_sample" | wc -l)

    # Filter out matching positions from the tumor VCF based on blood positions
    grep -vFf "$BLOOD_POSITIONS" "$tumor_sample" > "$output_file"

    # Pocet vyfiltrovanych datovych riadkov
    filtered_lines=$(grep -v "^#" "$output_file" | wc -l)

    # Calculate the number of filtered lines
    lines_deleted=$((total_lines - filtered_lines))

    # Output summary to CSV file
    echo "$(basename "$tumor_sample"),$total_lines,$filtered_lines,$lines_deleted" >> "$summary_file"

    # Cleanup temporary file
    rm "$BLOOD_POSITIONS"

    echo "Processed ${tumor_sample}: TotalLines=$total_lines, FilteredLines=$filtered_lines, LinesDeleted=$lines_deleted"
done

echo "Processing complete. Summary saved to $summary_file."


