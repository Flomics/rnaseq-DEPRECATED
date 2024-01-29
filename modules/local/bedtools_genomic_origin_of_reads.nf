process BEDTOOLS_GENOMIC_ORIGIN_OF_READS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    // def args = task.ext.args   ?: ''
    // prefix   = task.ext.prefix ?: "${meta.id}"
    // def paired_end = meta.single_end ? '' : '-pe'
    // def memory     = task.memory.toGiga() + "G"

    // def strandedness = 'non-strand-specific'
    // if (meta.strandedness == 'forward') {
    //     strandedness = 'strand-specific-forward'
    // } else if (meta.strandedness == 'reverse') {
    //     strandedness = 'strand-specific-reverse'
    // }
    '''
    awk '$3 == "exon"' !{gtf} > filtered_annotation_exon.gtf

    #bedtools and samtools to extract reads mapping to exons
    bedtools intersect -abam !{bam} -b filtered_annotaion_exon.gtf -wa -u | samtools view -F 0x04 - | cut -f 1 | sort | uniq -c > reads_in_exons.txt

    samtools view -L <(awk '$3 == "exon" {print $1 "\t" $4 "\t" $5}' filtered_annotation.gtf) -b sub_SRR8492578.markdup.sorted.bam > filtered_reads_no_exon.bam
    '''
    //unset DISPLAY
    // mkdir tmp
    // export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp



    // qualimap \\
    //     --java-mem-size=$memory \\
    //     rnaseq \\
    //     $args \\
    //     -bam $bam \\
    //     -gtf $gtf \\
    //     -p $strandedness \\
    //     $paired_end \\
    //     -outdir $prefix

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    // END_VERSIONS
}
