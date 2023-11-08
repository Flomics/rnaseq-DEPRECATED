process FLOMICS_TRACKHUBS{
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    tuple val(meta), path(bam)
    tuple val(meta), path(bai)

    output:
    path "trackhub/*_trackhub_links.tsv",     emit: trackhub_links
    path "trackhub/*_trackDb.txt",            emit: trackhub_trackDb_files
    path "assembly_hub/*_trackhub_links.tsv", emit: assembly_hub_links
    path "assembly_hub/*_trackDb.txt",        emit: assembly_hub_trackDb_files

    script:
    prefix  = task.ext.prefix ?: "${meta.id}"
    uuid = params.uuid
    project= params.project
    profile= workflow.profile

    """
    make_trackDB.sh $bam $prefix $project $uuid $profile
    cat trackhub/UCSC.txt >> trackhub/${prefix}_trackhub_links.tsv
    cat assembly_hub/UCSC.txt >> assembly_hub/${prefix}_trackhub_links.tsv
    """
}
