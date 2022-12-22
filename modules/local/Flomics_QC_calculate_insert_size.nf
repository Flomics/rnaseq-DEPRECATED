process FLOMICS_QC_CALCULATE_INSERT_SIZE{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    
    output:
    path("*.insert_size_median.tsv")   , emit: insert_size



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"


    """
    picard_transcriptome_insert.sh -b $bam -p $prefix
    """
}