#!/usr/bin/env python

import csv
import sys

if len(sys.argv) != 3:
    print("Usage: script.py input_file.tsv output_file.yml")
    sys.exit(1)

tsv_file = sys.argv[1]
yaml_file = sys.argv[2]

def convert_to_numeric(value):
    try:
        return float(value)
    except ValueError:
        return value

with open(tsv_file, mode='r', newline='') as file:
    reader = csv.DictReader(file, delimiter='\t')

    for row in reader:
        # Convert all values to numeric if possible
        data_row = {k: convert_to_numeric(v) for k, v in row.items()}
        break

output_str = "{"
output_str += ", ".join(f"{k}: {v}" for k, v in data_row.items())
output_str += "}"

with open(yaml_file, 'w') as f:
    f.write(output_str)
