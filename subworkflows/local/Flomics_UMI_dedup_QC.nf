//
// Flomics UMI-dedup QC
//

include { FLOMICS_QC_CALCULATE_UMI_DEDUP_RATE      } from '../../modules/local/Flomics_QC_UMI_dedup.nf'

workflow FLOMICS_UMI_DEDUP_QC{
    take:
    
    bam_files // channel: [ val(meta), [ bam_dedup ], [ bam_without_deduplication ]]

    main:

    ///
    /// Calculate the UMI deduplication rate
    ///
    FLOMICS_QC_CALCULATE_UMI_DEDUP_RATE(bam_files)
    ch_Flomics_UMI_dedup_rate_QC       = FLOMICS_QC_CALCULATE_UMI_DEDUP_RATE.out.umi_dedup_rate.collect()

    emit:
    umi_dedup_rate                       = FLOMICS_QC_CALCULATE_UMI_DEDUP_RATE.out.umi_dedup_rate                     // channel: [ val(meta), results_dir ]

}
