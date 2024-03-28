#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
table1= data.frame(matrix(ncol = 66, nrow = 0))
# names(table1)= c("Ala_tRNA_fc", "Arg_tRNA_fc", "Asn_tRNA_fc", "Asp_tRNA_fc", "Cys_tRNA_fc", "Gln_tRNA_fc", "Glu_tRNA_fc", "Gly_tRNA_fc", "His_tRNA_fc", "IG_C_gene_fc", "IG_C_pseudogene_fc",
#                      "IG_D_gene_fc", "IG_J_gene_fc", "IG_J_pseudogene_fc", "IG_V_gene_fc", "IG_V_pseudogene_fc", "IG_pseudogene_fc", "Ile_tRNA_fc", "Leu_tRNA_fc", "Lys_tRNA_fc", "Met_tRNA_fc",
#                      "Mt_rRNA_fc", "Mt_tRNA_fc", "Phe_tRNA_fc", "Pro_tRNA_fc", "Pseudo_tRNA_fc", "SeC(e)_tRNA_fc", "SeC_tRNA_fc", "Ser_tRNA_fc", "Sup_tRNA_fc", "TEC_fc", "TR_C_gene_fc", "TR_D_gene_fc",
#                      "TR_J_gene_fc", "TR_J_pseudogene_fc", "TR_V_gene_fc", "TR_V_pseudogene_fc", "Thr_tRNA_fc", "Trp_tRNA_fc", "Tyr_tRNA_fc", "Undet_tRNA_fc", "Val_tRNA_fc", "lncRNA_fc", "miRNA_fc",
#                      "misc_RNA_fc", "polymorphic_pseudogene_fc", "processed_pseudogene_fc", "protein_coding_fc", "pseudogene_fc", "rRNA_fc", "rRNA_pseudogene_fc", "ribozyme_fc", "sRNA_fc", "scRNA_fc",
#                      "scaRNA_fc", "snRNA_fc", "snoRNA_fc", "spike_in_fc", "transcribed_processed_pseudogene_fc", "transcribed_unitary_pseudogene_fc", "transcribed_unprocessed_pseudogene_fc",
#                      "translated_processed_pseudogene_fc", "translated_unprocessed_pseudogene_fc", "unitary_pseudogene_fc", "unprocessed_pseudogene_fc", "vault_RNA_fc")

table2= read.csv("multiqc_data/multiqc_featurecounts_biotype_plot.txt", header = T, row.names=1, sep= "\t")
names(table2) <- paste(names(table2), "_fc", sep="")

#table_final= bind_rows(table1,table2)
table_final <- table2
table_final[is.na(table_final)] <- 0
table_final= setDT(table_final, keep.rownames = "Sample")[]
write.table(table_final, "tmp.biotype_table.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
