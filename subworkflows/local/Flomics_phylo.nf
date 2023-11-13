//
// Flomics Phylogenetic classification
//

include { SEQKIT } from '../../modules/local/seqkit.nf'
include { KRAKEN2 } from '../../modules/local/kraken2.nf'
include { KRONA } from '../../modules/local/krona.nf'
include { UNTAR as UNTAR_KRAKEN2_DB } from '../../modules/nf-core/modules/untar/main.nf'


workflow FLOMICS_PHYLO{
    take:
    reads


    main:

    ch_reads = reads
    if (!params.phylo_skip_subsampling) {
        //
        // MODULE: Subsample the trimmed fastq files with SEQKIT
        //
        SEQKIT (
            ch_reads
        )
        ch_reads = SEQKIT.out.reads
    }

    //
    // MODULE: Untar kraken2_db
    //
    UNTAR_KRAKEN2_DB ( [ [:], params.kraken2_db ])
    ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { it[1] }

    //
    // MODULE: Perform kraken2
    //
    //ch_kraken2_db = Channel.fromPath(params.kraken2_db, checkIfExists: true)
    KRAKEN2 (
        ch_reads, ch_kraken2_db
    )
    KRAKEN2.out.report.map { meta, report -> [ report ] }.collect()

    //
    // MODULE: krona plot the kraken2 reports
    //
    KRONA (
        KRAKEN2.out.report.map { meta, report -> [ report ] }.collect()
    )

    emit:
    kraken2_report = KRAKEN2.out.report.map { meta, report -> [ report ] }.collect()
}
