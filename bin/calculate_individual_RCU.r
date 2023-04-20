#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(biomaRt)
# library(tools)

coverage_table= read.table(args[1], as.is = T)

names(coverage_table)= c("transcript", "base", "coverage")
coverage_table$transcript= gsub('\\.[^\\.]*$', '', coverage_table$transcript)

AUC_table_all_transcripts= data.frame("Transcript_id"= unique(coverage_table$transcript), "Read_coverage_uniformity_score"= "NA", "Gene_id"= "NA")

counter= 1
for(gene in unique(coverage_table$transcript)){
  unique_gene_coverage_table= coverage_table[coverage_table$transcript==gene,]
  AUC_transcript= sum(unique_gene_coverage_table$coverage)
  total_area= max(unique_gene_coverage_table$coverage)*length(unique_gene_coverage_table$coverage)
  Percentage_AUC_total= AUC_transcript/total_area*100
  AUC_table_all_transcripts[counter,2]= Percentage_AUC_total
  counter= counter+1
}

mart <- useMart("ensembl","hsapiens_gene_ensembl")
ensemble2gene <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                       filters = "ensembl_transcript_id",
                       values = AUC_table_all_transcripts$Transcript,
                       mart = mart)

counter= 1
for(transcript in AUC_table_all_transcripts$Transcript_id){
  gene_name = ensemble2gene[ensemble2gene$ensembl_transcript_id==transcript,3]
  if(length(gene_name>0)){
    AUC_table_all_transcripts$Gene_id[counter]= gene_name
  }
  counter= counter+1
}

AUC_table_all_transcripts= AUC_table_all_transcripts[order(AUC_table_all_transcripts$Read_coverage_uniformity_score, decreasing = T),]
write.table(AUC_table_all_transcripts, file= paste(args[2],"_gene_coverage_profile_table.tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
