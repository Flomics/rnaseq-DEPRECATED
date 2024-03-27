#!/bin/bash

# Check if less than 2 arguments are passed
if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <header_file> <output_file>"
  exit 1
fi

header_file="$1"

output_file="$2"

temp_file=$(mktemp)

cp "$header_file" "$temp_file"

# Loop through each *_mqc.yaml file
for mqc_file in *_mqc.yaml; do
    content=$(<"$mqc_file")
    echo "  $content" >> "$temp_file"
done

# Save to the output file
cp "$temp_file" "$output_file"