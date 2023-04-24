process FLOMICS_QC_CALCULATE_INDIVIDUAL_RCU{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    path transcript_to_gene_id_tsv
    
    output:
    path("*_gene_coverage_profile_table.tsv")   , emit: individual_RCUs



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    samtools sort $bam > ${prefix}_sorted.bam
    samtools depth -a ${prefix}_sorted.bam > ${prefix}_coverage.tsv
    calculate_individual_RCU.r ${prefix}_coverage.tsv $transcript_to_gene_id_tsv $prefix
    """
}