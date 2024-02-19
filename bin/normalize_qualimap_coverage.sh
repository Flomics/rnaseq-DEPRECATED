#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Input file
input_file="$1"

# Output file
output_file="$2"

# Print header to output file
head -n 1 "$input_file" > "$output_file"

# Calculate row sums and normalize each element
awk 'BEGIN {FS=OFS="\t"} NR>1 {
    sum=0;
    for(i=2; i<=NF; i++) {
        sum+=$i;
    }
    for(i=2; i<=NF; i++) {
        $i = $i / sum;
    }
    print;
}' "$input_file" >> "$output_file"
