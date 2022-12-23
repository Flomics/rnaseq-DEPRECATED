process FLOMICS_QC_PARSER{
    tag "$meta.id"
    label 'process_high'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    path multiqc_data

    output:
    path("*_fastqc_QC.tsv")   , emit: fastqc_files



    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def num_ends_sequencing = meta.single_end ? "single_end" : "paired_end"
    
    """
    fastQC_parser.sh $num_ends_sequencing $prefix #Parses the fastQC
    """
}