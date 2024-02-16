#!/bin/bash

# Define the header YAML file
header_file="genomic_origin_of_reads_header.yaml"

# Define the output file
output_file="genomic_origin_of_reads_mqc.yaml"

# Create an empty temporary file to store data
temp_file=$(mktemp)

# Copy header.yaml to the temporary file
cp "$header_file" "$temp_file"

# Loop through each *_mqc.yaml file
for mqc_file in *_mqc.yaml; do
    # Read the content of the mqc file
    content=$(<"$mqc_file")

    # Append the content to the temporary file
    echo "  $content" >> "$temp_file"
done

# Append the closing part to the temporary file
echo "}" >> "$temp_file"

# Save to the output file
cp "$temp_file" "$output_file"


