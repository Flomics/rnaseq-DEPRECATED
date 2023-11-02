process FLOMICS_QC_CALCULATE_LIBRARY_BALANCE{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    path(results_dir)   
    
    output:
    path("*_genes_contributing_to_percentage_reads.tsv")   , emit: library_balance_table



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"


    """
    gene_relative_abundance.r $prefix
    """
}