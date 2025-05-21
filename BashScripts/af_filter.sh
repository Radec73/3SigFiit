#!/bin/bash

# Input and output directories
input_folder="/home/cyprich/data/vqsr_filtered" # doplnit
output_folder="/home/cyprich/data/af_filtered" # doplnit
summary_file="/home/cyprich/data/af_filtered/summary.csv" # doplnit


# Create or overwrite the summary CSV file
echo "FileName,TotalLines,FilteredLines,LinesDeleted" > "$summary_file"

# Range for filtering AF
min_af=0.02
max_af=0.3

# Process each VCF file in the input directory
for input_file in "$input_folder"/*.vcf; do
    # Extract the file name without the path
    filename=$(basename "$input_file")
    output_file="$output_folder/$filename"

    # Total number of data lines (excluding header lines)
    total_lines=$(grep -v "^#" "$input_file" | wc -l)

    # Filter based on AF calculated from AD field
    awk -v min_af="$min_af" -v max_af="$max_af" '
    BEGIN { FS="\t"; OFS="\t" }
    {
        # Skip header lines and print them unchanged
        if ($0 ~ /^#/) {
            print $0
            next
        }

        # Extract the INFO field (10th column) and split the GT:AD:DP:GQ:PL field
        split($10, format_fields, ":")
        ad_field = format_fields[2]

        # Split AD field into reference and alternate allele depths
        split(ad_field, ad_values, ",")
        ref_depth = ad_values[1]
        alt_depth = ad_values[2]

        # Calculate allele frequency: AF = alt_depth / (ref_depth + alt_depth)
        if (ref_depth + alt_depth > 0) {
            af = alt_depth / (ref_depth + alt_depth)
        } else {
            af = 0  # Avoid division by zero if there is no depth for either allele
        }

        # Filter based on the AF
        if (af >= min_af && af <= max_af) {
            print $0
        }
    }' "$input_file" > "$output_file"

    # Count filtered lines
    filtered_lines=$(grep -v "^#" "$output_file" | wc -l)

    # Count deleted lines
    lines_deleted=$((total_lines - filtered_lines))

    # Append results to the summary file
    echo "$filename,$total_lines,$filtered_lines,$lines_deleted" >> "$summary_file"

    echo "Processed $filename: TotalLines=$total_lines, FilteredLines=$filtered_lines, LinesDeleted=$lines_deleted"
done

echo "Processing complete. Summary saved to $summary_file."

