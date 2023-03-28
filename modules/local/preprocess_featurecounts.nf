process PREPROCESS_FEATURECOUNTS{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(biotype_counts)

    output:
    tuple val(meta), path("*.biotype_counts.tsv"), emit: biotype_counts_processed

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed -e '1,11d' $biotype_counts > processed.tsv
    ( echo $prefix && cat processed.tsv ) > ${prefix}.biotype_counts.tsv

    """
}
