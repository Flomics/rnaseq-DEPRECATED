process FLOMICS_QC_SPLICED_READS{
    tag "$meta.id"
    label 'process_high'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    path("*.splicedReads.stats.tsv")   , emit: splicedReads_QC
    path("*.spliceJunctions.stats.tsv")   , emit: spliceJunctions_QC



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"


    """
    Flomics_QC_spliced_reads.sh \\
    -b $bam \\
    -g $gtf \\
    -p $prefix
    """
}