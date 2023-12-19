process FLOMICS_QC_CALCULATE_LIBRARY_BALANCE{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"

    input:
    path(counts)

    output:
    path("genes_contributing_to_percentage_reads.tsv")   , emit: library_balance_table

    script:
    """
    gene_library_diversity.r $counts
    """
}
