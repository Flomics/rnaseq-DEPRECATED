#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(stringr)
library(dplyr)

# load featurecounts results
featurecounts <- read.table(text = gsub(" ", "\t", readLines("/home/ctuni/Downloads/10K_CRC014.biotype_counts.tsv"))) #change this in pipeline with the path to bin
colnames(featurecounts) <- featurecounts[1,]
featurecounts <- featurecounts[-1,]
