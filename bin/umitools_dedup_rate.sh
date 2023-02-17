#!/bin/bash

lines_dedup_bam=$(samtools view -F 4 $1| cut -f1 | sort | uniq | wc -l)
lines_bam=$(samtools view -F 4 $2| cut -f1 | sort | uniq  | wc -l)

echo -e "Number of mapped unique molecules\tPercentage of unique molecules" >${3}_UMI_dedup.tsv
echo -en "$lines_dedup_bam\t" >>${3}_UMI_dedup.tsv
echo "scale=3; $lines_dedup_bam / $lines_bam *100" | bc >>${3}_UMI_dedup.tsv
