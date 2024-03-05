#!/usr/bin/env python

import csv
import yaml
import sys

# Check if the correct number of arguments are passed
if len(sys.argv) != 3:
    print("Usage: script.py input_file.tsv output_file.yml")
    sys.exit(1)

tsv_file = sys.argv[1]
yaml_file = sys.argv[2]

with open(tsv_file, mode='r', newline='') as file:
    # Read the TSV file
    reader = csv.DictReader(file, delimiter='\t')

    for row in reader:
        data_row = dict(row)
        break

# Convert the dictionary to a YAML string
yaml_str = yaml.dump(data_row, sort_keys=False)

with open(yaml_file, 'w') as f:
    f.write(yaml_str)

print(f"YAML file '{yaml_file}' has been created.")
