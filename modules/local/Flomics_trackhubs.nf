process FLOMICS_TRACKHUBS{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)

    output:
    path "*_trackhub_links.tsv", emit: trackhubs_path
    path "*_trackDb.txt", emit: trackDb_files


    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    outdir  = params.outdir
    runName = workflow.runName

    """
    make_trackDB.sh $bam $outdir $prefix $runName
    cat UCSC.txt >> ${prefix}_trackhub_links.tsv
    """
}