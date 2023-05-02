process FLOMICS_QC_CREATE_GENE_TO_TRANSCRIPT_TABLE{
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"

    input:
    path gtf

    output:
    path "transcript_to_gene_id.tsv", emit: transcript_to_gene_id_tsv

    shell:

    '''
    #The first part selects the transcript_id and gene names for all of the genes from GENCODE
    cut -f 9 !{gtf} | cut -f 2,4 -d ";" | grep "gene_name" | sed "s/ transcript_id //g" | sed "s/; gene_name/\t/" |
    sed 's/"//g' | sort -k1,1 | uniq > transcript_to_gene_id.tsv
    #The first part selects the transcript_id and gene names for all of the spike-ins. If there are no spike-ins, it will not take anything
    cut -f 9 !{gtf} | grep -v "##" | grep -v "gene_name" | cut -f 2,1 -d ";" | sed "s/; transcript_id/\t/g" | sed "s/gene_id//" |
    sed 's/"//g'|sort -k1,1 | uniq | awk ' BEGIN { OFS="\t" }{ print $2,$1}' >> transcript_to_gene_id.tsv
    '''

}