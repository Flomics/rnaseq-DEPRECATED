process KRONA {
    label 'process_high'

    container "flomicsbiotech/krona:latest"

    input:
    path(report)

    output:
    path("*html"), emit: krona

    script:
    """
    ktUpdateTaxonomy.sh
    ktImportTaxonomy -t 5 -m 3 -o multi-krona.html *.report.txt
    """
}
