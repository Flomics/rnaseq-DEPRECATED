//
// Flomics QC
//


include { FLOMICS_TRACKHUBS                             } from '../../modules/local/Flomics_trackhubs.nf'
include { FLOMICS_QC_SPLICED_READS                      } from '../../modules/local/Flomics_QC_spliced_reads.nf'
include { FLOMICS_QC_PARSER                             } from '../../modules/local/Flomics_FastQC_parser.nf'
include { FLOMICS_QC_CALCULATE_INSERT_SIZE              } from '../../modules/local/Flomics_QC_calculate_insert_size.nf'
include { FLOMICS_QC_CREATE_GENE_TO_TRANSCRIPT_TABLE    } from '../../modules/local/create_gene_to_transcript_table.nf'
include { FLOMICS_QC_CALCULATE_INDIVIDUAL_RCU           } from '../../modules/local/Flomics_QC_calculate_individual_RCU.nf'
include { FLOMICS_QC_CALCULATE_LIBRARY_BALANCE          } from '../../modules/local/Flomics_QC_calculate_library_balance.nf'
include { FLOMICS_QC_GINI_INDEX                         } from '../../modules/local/Flomics_QC_gini_index.nf'
include { FLOMICS_QC_SPIKE_INS                          } from '../../modules/local/Flomics_QC_spike_ins.nf'
include { FLOMICS_QC_AGGREGATOR                         } from '../../modules/local/Flomics_QC_agreggator.nf'
include { MAKE_FILES_PUBLIC                             } from '../../modules/local/make_files_public.nf'
include { FLOMICS_QC_KNIT                               } from '../../modules/local/Flomics_QC_knit.nf'


workflow FLOMICS_QC{
    take:
    samplesheet
    multiqc_data   // channel: multiqc_data/*
    multiqc_report // channel: multiqc_report.html
    bam_genome     // channel: [ val(meta), [ bam ]]
    bam_genome_indices // channel: [ val(meta), [ bam ]]
    bam_transcriptome // channel: [ val(meta), [ bam ]]
    gtf     // channel: /path/to/genome.gtf
    umi_dedup_rate_data
    salmon_results
    spike_in_concentration
    salmon_gene_tpm
    qc_dashboard
    salmon_gene_counts
    bedtools_genomic_origin_of_reads

    main:

    ///
    /// Makes the trackhubs and copies the bam files to the public bucket of s3
    ///
    ch_Flomics_trackhub_links         = Channel.empty()
    ch_Flomics_trackhub_trackDb_files = Channel.empty()
    ch_Flomics_assembly_hub_links     = Channel.empty()
    ch_Flomics_assembly_hub_trackDbs  = Channel.empty()
    FLOMICS_TRACKHUBS(bam_genome, bam_genome_indices)
    ch_Flomics_trackhub_links         = FLOMICS_TRACKHUBS.out.trackhub_links.collect()
    ch_Flomics_trackhub_trackDb_files = FLOMICS_TRACKHUBS.out.trackhub_trackDb_files.collect()
    ch_Flomics_assembly_hub_links     = FLOMICS_TRACKHUBS.out.assembly_hub_links.collect()
    ch_Flomics_assembly_hub_trackDbs  = FLOMICS_TRACKHUBS.out.assembly_hub_trackDb_files.collect()

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
    /// Create the gene_id to transcript tsv from the gtf file.
    ///
    ch_Flomics_gene_to_transcript       = Channel.empty()
    FLOMICS_QC_CREATE_GENE_TO_TRANSCRIPT_TABLE(gtf)
    ch_Flomics_gene_to_transcript       = FLOMICS_QC_CREATE_GENE_TO_TRANSCRIPT_TABLE.out.transcript_to_gene_id_tsv

    ///
    /// Calculate the RCU score for each gene individually
    ///
    ch_Flomics_individual_RCU_QC       = Channel.empty()
    FLOMICS_QC_CALCULATE_INDIVIDUAL_RCU(bam_transcriptome, ch_Flomics_gene_to_transcript)
    ch_Flomics_individual_RCU_QC       = FLOMICS_QC_CALCULATE_INDIVIDUAL_RCU.out.individual_RCUs.collect()

    ///
    /// Calculate the library balance (number of genes contributing to a percentage of the library)
    ///
    ch_Flomics_library_balance       = Channel.empty()
    FLOMICS_QC_CALCULATE_LIBRARY_BALANCE(salmon_gene_counts)
    ch_Flomics_library_balance       = FLOMICS_QC_CALCULATE_LIBRARY_BALANCE.out.library_balance_table

    ///
    /// Calculate gini index for the gene counts table
    ///
    ch_Flomics_gini_index       = Channel.empty()
    FLOMICS_QC_GINI_INDEX ( salmon_gene_counts )
    ch_Flomics_gini_index        = FLOMICS_QC_GINI_INDEX.out.gini_index_table.collect()

    ///
    /// Plot Spike-in concentration vs spike-in TPM and obtain correlation coefficients
    ///
    ch_Flomics_correlation_coefficients       = Channel.empty()
    FLOMICS_QC_SPIKE_INS ( spike_in_concentration, salmon_gene_tpm )
    ch_Flomics_correlation_coefficients        = FLOMICS_QC_SPIKE_INS.out.correlation_coefficients_table.collect()

    ///
    /// Aggregate all the QC from multiQC and extra QC into a new tsv
    ///
    ch_Flomics_QC_report            = Channel.empty()
    FLOMICS_QC_AGGREGATOR ( samplesheet, multiqc_data, ch_Flomics_trackhub_links, ch_Flomics_trackhub_trackDb_files, ch_Flomics_assembly_hub_links, ch_Flomics_assembly_hub_trackDbs, ch_Flomics_splicedReads_QC, ch_Flomics_spliceJunctions_QC,
    ch_Flomics_insert_size_QC, umi_dedup_rate_data, ch_Flomics_library_balance, ch_Flomics_FastQC, ch_Flomics_correlation_coefficients, bedtools_genomic_origin_of_reads)
    ch_Flomics_QC_report            = FLOMICS_QC_AGGREGATOR.out.flomics_report

    ///
    /// Publish the necessary files to the flomics-public bucket (same structure as the trackhubs runs)
    ///

    if (!workflow.profile.contains('test')){
        MAKE_FILES_PUBLIC (multiqc_report, ch_Flomics_QC_report)
    }
    ///
    /// Knit the Flomics QC into an interactive HTML dashboard
    ///
    FLOMICS_QC_KNIT ( FLOMICS_QC_AGGREGATOR.out.flomics_report, qc_dashboard )

}
