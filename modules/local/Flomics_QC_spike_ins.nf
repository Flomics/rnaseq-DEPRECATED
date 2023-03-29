process FLOMICS_QC_SPIKE_INS{
    tag "$meta"
    label 'process_low'

    container "flomicsbiotech/markdown_pkgs:dev"


    input:
    path(spike_in_concentration)
    path(gene_tpm)
    path(biotype_counts)

    output:
    path("correlation_coefs.tsv") , emit: correlation_coefficients_table
    path("*.png") , emit: correlation_png
    path("biotype_counts.tsv"), emit: biotype_counts_joined

    script:
    """
    plot_correlation.r $spike_in_concentration $gene_tpm
    paste $biotype_counts > biotype_counts.tsv

    """
    //spike-in_qc.r biotype_counts.tsv
}
