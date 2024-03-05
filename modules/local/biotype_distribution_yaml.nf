process BIOTYPE_DISTRIBUTION_YAML {
    label 'process_medium'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"

    input:
    path(yaml)
    path(header)

    output:
    path("biotype_distribution_mqc.yaml"), emit: mqc_biotype_distribution


    script:
    """
    mkdir tmp/
    cp $yaml tmp/
    cp $header tmp/
    cd tmp/
    bash make_yaml.sh $header $yaml
    cd ../
    cp tmp/biotype_distribution_mqc.yaml .
    """
}
