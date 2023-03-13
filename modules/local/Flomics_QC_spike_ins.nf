process FLOMICS_QC_SPIKE_INS{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(spike_in_concentration)
    path(gene_tpm)

    output:
    path("correlation_coefs.tsv") , emit: correlation_coefficients_table
    path("plots_correlation.pdf") , emit: correlation_pdf

    script:
    """
    plot_correlation.r $spike_in_concentration $gene_tpm
    """
}
