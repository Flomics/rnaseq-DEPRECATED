#!/bin/bash

# Check if less than 2 arguments are passed
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <header_file> <output_file>"
  exit 1
fi

# Define the header YAML file from the first command line argument
header_file="$1"

# Define the output file from the second command line argument
output_file="$2"

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
#echo "}" >> "$temp_file"

# Save to the output file
cp "$temp_file" "$output_file"

