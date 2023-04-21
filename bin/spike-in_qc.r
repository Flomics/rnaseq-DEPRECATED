#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(stringr)
library(dplyr)

# load featurecounts results
featurecounts <- read.table(text = gsub(" ", "\t", readLines("/home/ctuni/Downloads/biotype_counts(3).tsv"))) #change this in pipeline with the path to bin
colnames(featurecounts) <- featurecounts[1,]
featurecounts <- featurecounts[-1,]
rownames(featurecounts) <- featurecounts[,1]

featurecounts <- featurecounts[ , -which(names(featurecounts) %in% c("biotype"))]
featurecounts_num <- mutate_all(featurecounts, function(x) as.numeric(x))

featurecounts_t <- t(featurecounts_num)

featurecounts_ordered <- featurecounts_t[order(rownames(featurecounts_t)),]
