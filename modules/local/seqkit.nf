process SEQKIT {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::seqkit=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/seqkit:2.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("sub_*"), emit: reads

    script:
    if (meta.single_end) {
        """
        zcat $reads | seqkit sample -p $params.phylo_proportion -s 80 -o sub_$reads
        """
    } else {
        """
        zcat ${reads[0]} | seqkit sample -p $params.phylo_proportion -s 80 -o sub_${reads[0]}
        zcat ${reads[1]} | seqkit sample -p $params.phylo_proportion -s 80 -o sub_${reads[1]}
        """
    }
}
