#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

library(edgeR)

seqdata <- read.delim(args[1], stringsAsFactors = FALSE)
countdata <- seqdata[, -(1:2)]

g <- as.data.frame(gini(countdata))

g <- data.frame(rownames(g), g)
colnames(g) <- c("samplename", "gini_index")

write.table(g, file = "gini_index.tsv", row.names = FALSE, col.names = TRUE)
