#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
salmon_quant= read.csv(paste(args[1],"/quant.sf", sep = ""), header = T, sep= "\t")

salmon_quant= salmon_quant[order(salmon_quant$TPM, decreasing = T),]

genes_1_percent= length(which(cumsum(salmon_quant$TPM) <= 10000))+1
genes_5_percent= length(which(cumsum(salmon_quant$TPM) <= 50000))+1
genes_10_percent= length(which(cumsum(salmon_quant$TPM) <= 100000))+1
genes_50_percent= length(which(cumsum(salmon_quant$TPM) <= 500000))+1
genes_80_percent= length(which(cumsum(salmon_quant$TPM) <= 800000))+1

gene_calculation_table= data.frame("Genes contributing to 1% of reads"= genes_1_percent, "Genes contributing to 5% of reads"= genes_5_percent,
                                   "Genes contributing to 10% of reads"= genes_10_percent,"Genes contributing to 50% of reads"= genes_50_percent,
                                   "Genes contributing to 80% of reads"= genes_80_percent, check.names = FALSE)

write.table(gene_calculation_table, file= paste(args[1],"_genes_contributing_to_percentage_reads.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = F)
