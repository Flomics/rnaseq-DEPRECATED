#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <bam_file> <gtf_file> <sample>"
    exit 1
fi

# Assign command line arguments to variables
bam=$1
gtf=$2
sample=$3

#take only exon records
awk '$3 == "exon"' !{gtf} > filtered_annotation_exon.gtf

#name sort input bam
samtools sort -n !{bam} -o name_sorted.bam

#extract gene locus coordinates
cat !{gtf} | extract_locus_coords.pl - > gencode.genes.bed

#mapped fragments
samtools view -F4 name_sorted.bam | cut -f1 | sort | uniq > "${sample}_mapped_fragments.list.txt"

#exonic fragments.
#Do not take strandedness into account for simplicity. If we did, we would need to take the library type into account (e.g. only read#2 should match the strand in SMARTEr libraries)
#Use -split to avoid taking into account gaps in spliced alignments
# takes about 1.3Gb of mem without the -sorted option
bedtools intersect -split -abam name_sorted.bam -b filtered_annotation_exon.gtf -wa -u  | samtools view -F4 - | cut -f1 | sort | uniq > "${sample}_exonic_fragments.list.txt"

#genic fragments
bedtools intersect -split -abam name_sorted.bam -b gencode.genes.bed -wa -u  | samtools view -F4 - | cut -f1 | sort | uniq > "${sample}_genic_fragments.list.txt"

#intergenic fragments (mapped reads that do not overlap genic regions)
comm -3 mapped_fragments.list.txt genic_fragments.list.txt > "${sample}_intergenic_fragments.list.txt"

#intronic fragments (genic fragments that do not overlap exons)
comm -3 genic_fragments.list.txt exonic_fragments.list.txt  > "${sample}_intronic_fragments.list.txt"
