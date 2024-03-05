process BIOTYPE_DISTRIBUTION_YAML {
    label 'process_medium'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"

    input:
    path(tsv)
    path(header)

    output:
    path("biotype_distribution_mqc.yaml"), emit: mqc_biotype_distribution


    script:
    """
    tsv_to_yaml.py $tsv tmp.yaml
    mkdir tmp/
    cp tmp.yaml tmp/
    cp $header tmp/
    cd tmp/
    bash make_yaml.sh $header biotype_distribution_mqc.yaml
    cd ../
    cp tmp/biotype_distribution_mqc.yaml .
    """
}
