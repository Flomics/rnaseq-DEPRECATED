process UMITOOLS_EXTRACT {
    tag "$meta.id"
    label "process_low"

    container "flomicsbiotech/umitools_last_version:latest"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads
    tuple val(meta), path("*.log")     , emit: log
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        export PYTHONHASHSEED=0
        umi_tools \\
            extract \\
            --random-seed=123456789\\
            -I $reads \\
            -S ${prefix}.umi_extract.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
        END_VERSIONS
        """
    }  else {
        """
        export PYTHONHASHSEED=0
        umi_tools \\
            extract \\
            --random-seed=123456789\\
            -I ${reads[0]} \\
            --read2-in=${reads[1]} \\
            -S ${prefix}.umi_extract_1.fastq.gz \\
            --read2-out=${prefix}.umi_extract_2.fastq.gz \\
            $args \\
            > ${prefix}.umi_extract.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            umitools: \$(umi_tools --version 2>&1 | sed 's/^.*UMI-tools version://; s/ *\$//')
        END_VERSIONS
        """
    }
}
