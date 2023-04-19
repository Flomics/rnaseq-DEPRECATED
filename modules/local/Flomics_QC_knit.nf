process FLOMICS_QC_KNIT{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(flomics_report)
    path qc_dashboard


    output:
    path("QC_dashboard.html") , emit: qc_dashboard

    script:
    outdir  = params.outdir
    """
    cp $qc_dashboard flomics_qc_dashboard.Rmd
    Rscript -e 'flomics_report="${flomics_report}"; rmarkdown::render(input = "flomics_qc_dashboard.Rmd", output_file = "QC_dashboard.html")'
    """
}
