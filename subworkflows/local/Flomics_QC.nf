//
// Flomics QC
//


include { FLOMICS_TRACKHUBS                     } from '../../modules/local/Flomics_trackhubs.nf'
include { FLOMICS_QC_SPLICED_READS              } from '../../modules/local/Flomics_QC_spliced_reads.nf'
include { FLOMICS_QC_PARSER                     } from '../../modules/local/Flomics_FastQC_parser.nf'
include { FLOMICS_QC_CALCULATE_INSERT_SIZE      } from '../../modules/local/Flomics_QC_calculate_insert_size.nf'
include { FLOMICS_QC_CALCULATE_LIBRARY_BALANCE  } from '../../modules/local/Flomics_QC_calculate_library_balance.nf'
include { FLOMICS_QC_AGGREGATOR                 } from '../../modules/local/Flomics_QC_agreggator.nf'

workflow FLOMICS_QC{
    take:
    multiqc_data    // channel: multiqc_data/* 
    bam_genome     // channel: [ val(meta), [ bam ]]
    bam_genome_indices // channel: [ val(meta), [ bam ]]
    bam_transcriptome // channel: [ val(meta), [ bam ]]
    gtf     // channel: /path/to/genome.gtf
    umi_dedup_rate_data
    salmon_results


    main:

    ///
    /// Makes the trackhubs and copies the bam files to the public bucket of s3
    ///
    ch_Flomics_trackhubs            = Channel.empty()
    ch_Flomics_trackDbs            = Channel.empty()
    FLOMICS_TRACKHUBS(bam_genome, bam_genome_indices)
    ch_Flomics_trackhubs            = FLOMICS_TRACKHUBS.out.trackhubs_path.collect()
    ch_Flomics_trackDbs            = FLOMICS_TRACKHUBS.out.trackDb_files.collect()

    ///
    /// Calculate the percentage of spliced reads and the percentage of splice junctions
    ///
    ch_Flomics_splicedReads_QC      = Channel.empty()
    ch_Flomics_spliceJunctions_QC   = Channel.empty()
    FLOMICS_QC_SPLICED_READS( bam_genome, gtf)
    ch_Flomics_splicedReads_QC      =   FLOMICS_QC_SPLICED_READS.out.splicedReads_QC.collect()
    ch_Flomics_spliceJunctions_QC   =   FLOMICS_QC_SPLICED_READS.out.spliceJunctions_QC.collect()

    ///
    /// Parse the FastQC file
    ///
    ch_Flomics_FastQC       = Channel.empty()
    FLOMICS_QC_PARSER(bam_genome_indices, multiqc_data)
    ch_Flomics_FastQC       = FLOMICS_QC_PARSER.out.fastqc_files.collect()

    ///
    /// Calculate the insert size
    ///
    ch_Flomics_insert_size_QC       = Channel.empty()
    FLOMICS_QC_CALCULATE_INSERT_SIZE(bam_transcriptome)
    ch_Flomics_insert_size_QC       = FLOMICS_QC_CALCULATE_INSERT_SIZE.out.insert_size.collect()

    ///
    /// Calculate the library balance (number of genes contributing to a percentage of the library)
    ///
    ch_Flomics_library_balance       = Channel.empty()
    FLOMICS_QC_CALCULATE_LIBRARY_BALANCE(bam_transcriptome, salmon_results)
    ch_Flomics_library_balance       = FLOMICS_QC_CALCULATE_LIBRARY_BALANCE.out.library_balance_table.collect()

    ///
    /// Aggregate all the QC from multiQC and extra QC into a new tsv
    ///
    FLOMICS_QC_AGGREGATOR ( multiqc_data, ch_Flomics_trackhubs, ch_Flomics_trackDbs, ch_Flomics_splicedReads_QC, ch_Flomics_spliceJunctions_QC,
    ch_Flomics_insert_size_QC, umi_dedup_rate_data, ch_Flomics_library_balance, ch_Flomics_FastQC)

    
}