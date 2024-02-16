process GENOMIC_ORIGIN_OF_READS_YAML {
    tag "$meta"
    label 'process_medium'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"

    input:
    tuple val(meta), path(yaml)
    path(header)

    output:
    path("genomic_origin_of_reads_mqc.yaml"), emit: mqc_genomic_origin_of_reads


    script:
    """
    mkdir tmp/
    cp $yaml tmp/
    cp $header tmp/
    cd tmp/
    bash make_yaml.sh
    """
}
