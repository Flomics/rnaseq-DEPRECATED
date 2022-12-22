#!/bin/bash
while getopts "b:p:" arg; do
  case $arg in
    b) bam=$OPTARG;;
    p) prefix=$OPTARG;;

  esac
done

: ${bam?"No bam file provided through -b"}
: ${prefix?"No prefix provided through -p"}

samtools sort $bam > ${prefix}.Aligned.toTranscriptome.sorted.bam
picard CollectInsertSizeMetrics I=${prefix}.Aligned.toTranscriptome.sorted.bam O=${prefix}_picard_metrics.txt H=${prefix}.insert_size_histogram.pdf M=0.5
cat ${prefix}_picard_metrics.txt | grep -A1 "MEDIAN" | cut -f1 > ${prefix}.insert_size_median.tsv
