#!/usr/bin/env Rscript

library(pROC)

gene_coverage_matrix= read.csv("multiqc_data/mqc_qualimap_gene_coverage_profile_Normalised.txt", sep = "\t", row.names = 1)
gene_coverage_profile= data.frame("Sample"= rownames(gene_coverage_matrix), "gene_coverage_profile"= 1:length(rownames(gene_coverage_matrix)))

for (n in  1:length(rownames(gene_coverage_matrix))){
  gene_coverage_sample= gene_coverage_matrix[n,]
  df= data.frame("gene_coverage_sample"= as.numeric(gene_coverage_sample), "gene_coverage_point"= c(0:99))
  AUC_sample= auc(df$gene_coverage_point, df$gene_coverage_sample)
  total_area= max(df$gene_coverage_sample)*100
  Percentage_AUC_total= AUC_sample/total_area*100
  gene_coverage_profile$gene_coverage_profile[n]= Percentage_AUC_total
}

write.table(gene_coverage_profile, file= "gene_coverage_profile_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
