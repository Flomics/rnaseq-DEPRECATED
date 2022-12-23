process FLOMICS_QC_CALCULATE_UMI_DEDUP_RATE{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bam_dedup)
    
    output:
    path("*_UMI_dedup.tsv")   , emit: umi_dedup_rate



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"


    """
    umitools_dedup_rate.sh $bam_dedup $bam $prefix
    """
}