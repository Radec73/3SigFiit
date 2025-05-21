#!/bin/bash

# Define paths
input_folder="/home/cyprich/data/filtered_Artifacts"
output_folder="/home/cyprich/data/VQSR_processed"
summary_file="/home/cyprich/data/summary_files/vqsr_summary.csv"
vqsr_output_folder="/home/cyprich/data/VQSR_output/"
gatk_path="/opt/gatk-4.4.0.0/gatk"
resources="/home/cyprich/data/resources/"

echo "Input File,Total Lines,Passed Lines,To be Deleted Lines" > "$summary_file"

# Loop through all .vcf files in the input directory
for input_file in ${input_folder}/*.vcf; do
    # Extract the base file name without directory and extension
    file_name=$(basename "$input_file" .vcf)
    # Define output file names
    output_file="${output_folder}/${file_name}_vqsr.vcf"
    # Get the total number of data lines (non-header) in the input file
    total_lines=$(grep -v "^#" "$input_file" | wc -l)
    # Run VariantRecalibrator
    python3 ${gatk_path} VariantRecalibrator \
       -R ${resources}hg38.fa \
       -V ${input_file} \
       --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${resources}hapmap_3.3.hg38.vcf.gz \
       --resource:omni,known=false,training=true,truth=false,prior=12.0 ${resources}1000G_omni2.5.hg38.vcf.gz \
       --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${resources}1000G_phase1.snps.high_confidence.hg38.vcf.gz \
       --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${resources}Homo_sapiens_assembly38.dbsnp138.vcf \
       -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
       --mode SNP \
       -O ${vqsr_output_folder}output.recal \
       --tranches-file ${vqsr_output_folder}output.tranches

    # Run ApplyVQSR
    python3 ${gatk_path} ApplyVQSR \
       -R ${resources}hg38.fa \
       -V ${input_file} \
       -O ${output_file} \
       -ts-filter-level 90 \
       --tranches-file ${vqsr_output_folder}output.tranches \
       --recal-file ${vqsr_output_folder}output.recal \
       --mode SNP
       
    passed_lines=$(awk 'BEGIN {FS="\t"} !/^#/ && $7 == "PASS" {count++} END {print count}' "$output_file")

   # Calculate the number of lines deleted (total lines - PASS lines)
    not_passed_lines=$((total_lines - passed_lines))

    # Output summary to CSV file
    echo "$file_name,$total_lines,$filtered_lines,$not_passed_lines" >> $summary_file

    echo "Processed $file_name"
done

echo "Processing complete. Summary saved to $summary_file."
