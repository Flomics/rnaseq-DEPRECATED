process FLOMICS_QC_GINI_INDEX{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(gene_counts)

    output:
    path("gini_index.tsv") , emit: gini_index_table

    script:
    """
    gini_index.r $gene_counts

    """
}
