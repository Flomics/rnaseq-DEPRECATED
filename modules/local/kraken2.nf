process KRAKEN2 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::kraken2=2.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_3' :
        'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_3' }"

    input:
    tuple val(meta), path(reads)
    path(db)

    output:
    tuple val(meta), path("*report.txt"), emit: report

    script:
    pe = meta.single_end ? "" : "--paired"
    if (meta.single_end) {
        """
        kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${meta.id}.kraken2.report.txt \\
        $pe \\
        --gzip-compressed \\
        $reads
        """
    } else {
        """
        kraken2 \\
        --db $db \\
        --threads $task.cpus \\
        --report ${meta.id}.kraken2.report.txt \\
        $pe \\
        --gzip-compressed \\
        ${reads[0]} ${reads[1]}
        """
    }

}
