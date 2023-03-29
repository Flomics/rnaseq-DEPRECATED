process MULTIQC_CUSTOM_BIOTYPE {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(count)
    path  header

    output:
    tuple val(meta), path("*biotype_counts_mqc.tsv"), emit: biotype_counts
    path "versions.yml"           , emit: versions
    tuple val(meta), path("*.tsv"), emit: tsv

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cut -f 1,7 $count | tail -n +3 | cat $header - >> ${prefix}.biotype_counts_mqc.tsv

    mqc_features_stat.py \\
        ${prefix}.biotype_counts_mqc.tsv \\
        -s $meta.id \\
        -f rRNA \\
        -o ${prefix}.biotype_counts_rrna_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
