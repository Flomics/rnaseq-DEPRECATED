#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
table1= data.frame(matrix(ncol = 66, nrow = 0))
names(table1)= c("Ala_tRNA", "Arg_tRNA", "Asn_tRNA", "Asp_tRNA", "Cys_tRNA", "Gln_tRNA", "Glu_tRNA", "Gly_tRNA", "His_tRNA", "IG_C_gene", "IG_C_pseudogene",
                     "IG_D_gene", "IG_J_gene", "IG_J_pseudogene", "IG_V_gene", "IG_V_pseudogene", "IG_pseudogene", "Ile_tRNA", "Leu_tRNA", "Lys_tRNA", "Met_tRNA",
                     "Mt_rRNA", "Mt_tRNA", "Phe_tRNA", "Pro_tRNA", "Pseudo_tRNA", "SeC(e)_tRNA", "SeC_tRNA", "Ser_tRNA", "Sup_tRNA", "TEC", "TR_C_gene", "TR_D_gene",
                     "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene", "Thr_tRNA", "Trp_tRNA", "Tyr_tRNA", "Undet_tRNA", "Val_tRNA", "lncRNA", "miRNA",
                     "misc_RNA", "polymorphic_pseudogene", "processed_pseudogene", "protein_coding", "pseudogene", "rRNA", "rRNA_pseudogene", "ribozyme", "sRNA", "scRNA",
                     "scaRNA", "snRNA", "snoRNA", "spike_in", "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene", "transcribed_unprocessed_pseudogene",
                     "translated_processed_pseudogene", "translated_unprocessed_pseudogene", "unitary_pseudogene", "unprocessed_pseudogene", "vault_RNA")
table2= read.csv("multiqc_data/mqc_featurecounts_biotype_plot.txt", header = T, row.names=1, sep= "\t")

table_final= bind_rows(table1,table2)
table_final[is.na(table_final)] <- 0
table_final= setDT(table_final, keep.rownames = "Sample")[]
write.table(table_final, "biotype_table.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
