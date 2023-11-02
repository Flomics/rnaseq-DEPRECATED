process FLOMICS_QC_CALCULATE_INSERT_SIZE{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    
    output:
    path("*.insert_size_median.tsv")   , emit: insert_size
    path("*.insert_size_histogram.pdf"), optional:true, emit: histogram_pdf


    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    def num_ends_sequencing = meta.single_end ? "single_end" : "paired_end"

    """
    if [[ "$num_ends_sequencing" == "paired_end" ]]
    then
        picard_transcriptome_insert.sh -b $bam -p $prefix
    else
        echo -e "MEDIAN_INSERT_SIZE\nNA" > ${prefix}.insert_size_median.tsv
    fi
    """
}
