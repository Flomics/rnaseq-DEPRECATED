process CONCAT_FEATURECOUNTS{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(biotype_counts)

    output:
    path("biotype_counts.tsv"), emit: biotype_counts_joined

    script:
    """
    paste $biotype_counts > biotype_counts.tsv
    """
}
