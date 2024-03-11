#!/usr/bin/env python3

import pandas as pd
import sys

def transpose_tsv(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')

    df_transposed = df.T

    df_transposed.to_csv(output_file, sep='\t', header=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file.tsv output_file.tsv")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]

        transpose_tsv(input_file, output_file)
