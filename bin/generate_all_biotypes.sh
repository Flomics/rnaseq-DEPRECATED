#!/bin/bash

while getopts "g:" arg; do
  case $arg in
    g) gtf=$OPTARG;;
  esac
done

: ${gtf?"No gtf file provided through -g"}

cat $gtf | grep "gene_type\|gene_biotype" |   awk -F'gene_type|gene_biotype' '{ print $2 "  " }' | cut -d'"' -f2 | sort | uniq > gene_type.txt
tr "\n" "\t" < gene_type.txt > gene_type.tsv
sed -i 's/[ \t]*$//' gene_type.tsv
