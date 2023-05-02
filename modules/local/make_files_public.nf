process MAKE_FILES_PUBLIC{
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    path multiqc_report
    path flomics_QC_report_table


    output:

    script:
    timestamp = workflow.start
    project= params.project


    """
    make_files_public.sh $multiqc_report $flomics_QC_report_table $timestamp $project

    """
}