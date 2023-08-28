process FLOMICS_QC_GINI_INDEX{
    tag "$meta"
    label 'process_low'

    conda "bioconda::bioconductor-edger=3.40.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-edger:3.40.0--r42hf17093f_1' :
        'quay.io/biocontainers/bioconductor-edger:3.40.0--r42hf17093f_1' }"



    input:
    path(gene_counts)

    output:
    path("gini_index.tsv") , emit: gini_index_table

    script:
    """
    gini_index.r $gene_counts
    """
}
