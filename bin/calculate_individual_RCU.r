#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

coverage_table= read.table(args[1], as.is = T)
transcript_to_gene_table= read.table(args[2], as.is = T, header = T)

names(coverage_table)= c("transcript", "base", "coverage")

AUC_table_all_transcripts= data.frame("Transcript_id"= unique(coverage_table$transcript), "Read_coverage_uniformity_score"= "NA", "Gene_id"= "NA")

#Index the transcripts (length of the transcripts)
table_index_genes= table(coverage_table$transcript)
sorted_table_index_genes = table_index_genes[unique(coverage_table$transcript)]
indices= list()
start_num=1
for (i in 1:length(sorted_table_index_genes)){
  indices[[i]]= c(start_num,start_num+sorted_table_index_genes[[i]]-1)
  start_num= start_num+sorted_table_index_genes[[i]]
}

#Calculate the RCU for each transcript and select the gene id that corresponds to the transcript
pre_time= Sys.time()
for(gene_num in 1:length(unique(coverage_table$transcript))){
  unique_gene_coverage_table= coverage_table[c(indices[[gene_num]][1]:indices[[gene_num]][2]),]
  AUC_transcript= sum(unique_gene_coverage_table$coverage)
  total_area= max(unique_gene_coverage_table$coverage)*length(unique_gene_coverage_table$coverage)
  Percentage_AUC_total= AUC_transcript/total_area*100
  AUC_table_all_transcripts$Read_coverage_uniformity_score[gene_num]= Percentage_AUC_total
  gene_name = transcript_to_gene_table[transcript_to_gene_table$Transcript_id==unique(unique_gene_coverage_table$transcript),2]
  AUC_table_all_transcripts$Gene_id[gene_num]= gene_name
}
run_time= Sys.time() - pre_time
print(run_time)

AUC_table_all_transcripts= AUC_table_all_transcripts[order(as.numeric(AUC_table_all_transcripts$Read_coverage_uniformity_score), decreasing = T),]
write.table(AUC_table_all_transcripts, file= paste(args[3],"_gene_coverage_profile_table.tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
