process FLOMICS_QC_AGGREGATOR{
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"


    input:
    path samplesheet
    path multiqc_data
    path "trackhub/*"
    path "trackhub/*"
    path "assembly_hub/*"
    path "assembly_hub/*"
    path splicedReads_QC
    path spliceJunctions_QC
    path insert_size
    path umi_dedup_rate_data
    path library_balance_data
    path fastqc_files
    path spike_in_qc

    output:
    path "QC_table.tsv", emit: flomics_report


    shell:

    project = params.project
    uuid = params.uuid
    profile = workflow.profile

    '''

    cut -f1 -d "," !{samplesheet} | sed -e "s/sample/Sample/" | head -n 1 > samplenames.tsv
    cut -f1 -d "," !{samplesheet} | tail -n+2 | sort | uniq >> samplenames.tsv

    tail -n +2 samplenames.tsv > samples.tsv
    merge_trackDB.sh !{project} !{uuid} !{profile} #Merges the trackDb files and uploads it to s3

    echo -e "Sample\ttrackhub_link" > trackhub/trackhub_links.tsv
    while IFS= read -r sample; do
        if test -f "trackhub/${sample}_trackhub_links.tsv"; then
            cat trackhub/${sample}_trackhub_links.tsv | sed "s/^/$sample\t/" >> trackhub/trackhub_links.tsv
        else
            echo -e "$sample\tNA" >> trackhub/trackhub_links.tsv
        fi
    done < samples.tsv

    echo -e "Sample\ttrackhub_link" > assembly_hub/trackhub_links.tsv
    while IFS= read -r sample; do
        if test -f "assembly_hub/${sample}_trackhub_links.tsv"; then
            cat assembly_hub/${sample}_trackhub_links.tsv | sed "s/^/$sample\t/" >> assembly_hub/trackhub_links.tsv
        else
            echo -e "$sample\tNA" >> assembly_hub/trackhub_links.tsv
        fi
    done < samples.tsv

    rm -f tmp.biotype_table.tsv
    biotype_table_parser.r  #Sorts the columns of the biotypes according to last table
    cat tmp.biotype_table.tsv | head -n 1 > biotype_table.tsv
    cat tmp.biotype_table.tsv | tail -n+2 | sort -k1,1 >> biotype_table.tsv

    #Normalize the qualimap provided coverage
    normalize_qualimap_coverage.sh multiqc_data/qualimap_rnaseq_cov_hist.txt multiqc_data/mqc_qualimap_gene_coverage_profile_Normalised.txt

    rm -f tmp.read_coverage_uniformity_score.tsv
    gene_coverage_profile_calculation.r #Calculates the gene coverage profile
    cat tmp.read_coverage_uniformity_score.tsv | head -n 1 > read_coverage_uniformity_score.tsv
    cat tmp.read_coverage_uniformity_score.tsv | tail -n+2 | sort -k1,1 >> read_coverage_uniformity_score.tsv


    echo -e "Sample\tper_base_sequence_quality\tper_sequence_quality_scores\tper_base_sequence_content\tper_sequence_gc_content\tper_base_n_content\tsequence_length_distribution\tsequence_duplication_levels\toverrepresented_sequences\tadapter_content" > fastqc_QC.tsv
    while IFS= read -r sample; do
        if test -f "${sample}_fastqc_QC.tsv"; then
            cat "${sample}_fastqc_QC.tsv" | sed "s/^/$sample\t/" >> fastqc_QC.tsv
        else
            echo -e "$sample\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA" >> fastqc_QC.tsv
        fi
    done < samples.tsv

    echo -e "Sample\tuniquely_mapped_reads\tspliced_reads\tpercentage_of_spliced_reads" > splicedReads_grouped.stats.tsv
    while IFS= read -r sample; do
        if test -f "${sample}.splicedReads.stats.tsv"; then
            tail -n +2 ${sample}.splicedReads.stats.tsv | sed "s/^/$sample\t/" >> splicedReads_grouped.stats.tsv
        else
            echo -e "$sample\tNA\tNA\tNA" >> splicedReads_grouped.stats.tsv
        fi
    done < samples.tsv

    echo -e "Sample\ttotal_splice_junctions\tknown_splice_junctions\t%known_splice_junctions" > spliceJunctions_grouped.stats.tsv
    while IFS= read -r sample; do
        if test -f "${sample}.spliceJunctions.stats.tsv"; then
            tail -n +2 ${sample}.spliceJunctions.stats.tsv | sed "s/^/$sample\t/" >> spliceJunctions_grouped.stats.tsv
        else
            echo -e "$sample\tNA\tNA\tNA" >> spliceJunctions_grouped.stats.tsv
        fi
    done < samples.tsv

    echo -e "Sample\tmedian_insert_size" > insert_size.tsv
    while IFS= read -r sample; do
        if test -f "${sample}.insert_size_median.tsv"; then
            tail -n +2 ${sample}.insert_size_median.tsv | sed "s/^/$sample\t/" >> insert_size.tsv
        else
            echo -e "$sample\tNA" >> insert_size.tsv
        fi
    done < samples.tsv

    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' genes_contributing_to_percentage_reads.tsv | head -n 1 > library_balance.tsv
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' genes_contributing_to_percentage_reads.tsv | tail -n+2 | sort -k1,1 >> library_balance.tsv

    echo -e "Sample\tNumber_mapped_unique_molecules\tPercentage_unique_molecules" > UMI_dedup_grouped.tsv
    while IFS= read -r sample; do
        if test -f "${sample}_UMI_dedup.tsv"; then
            tail -n +2 ${sample}_UMI_dedup.tsv| sed "s/^/$sample\t/" >> UMI_dedup_grouped.tsv
        else
            echo -e "$sample\tNA\tNA" >> UMI_dedup_grouped.tsv
        fi
    done < samples.tsv

    echo -e "Sample\tRead_number\tReads_passing_trimming\tPercentage_reads_passing_trimming"> cutadapt_QC.tsv
    cut -f1 multiqc_data/multiqc_cutadapt.txt | sed 's/_[0-9]$//g' > tmp.cutadapt.txt
    paste tmp.cutadapt.txt multiqc_data/multiqc_cutadapt.txt | awk '!seen[$1]++' > tmp2.cutadapt.txt
    cut -f 1,4,6 tmp2.cutadapt.txt | tail -n +2 | awk '{print $1"\t"$2"\t"$3"\t"($3/$2*100)}' >> tmp.cutadapt_QC.tsv
    cat tmp.cutadapt_QC.tsv | sort -k1,1 >> cutadapt_QC.tsv

    echo -e "Sample\tExonic\tIntronic\tIntergenic\tExonic_percentage\tIntronic_percentage\tIntergenic_percentage" > qualimap_QC.tsv
    tail -n +2 multiqc_data/qualimap_rnaseq_genome_results.txt | awk '{print $0"\t"$2/($2+$3+$4)*100"\t"$3/($2+$3+$4)*100"\t"$4/($2+$3+$4)*100'} | sort -k1,1 >> qualimap_QC.tsv

    echo -e "Sample\ttotal_reads\tavg_input_read_length\tnumber_of_uniquely_mapped_reads\tpercentage_of_uniquely_mapped_reads\tavg_mapped_read_length\tnumber_of_multimapped_reads\tpercentage_of_unmapped_too_short_reads\tmapped_percentage\taverage_mapped_length_percentage" > STAR_QC.tsv
    tail -n +2  multiqc_data/multiqc_star.txt | cut -f1,2,3,4,5,6,18,23 | awk '{print $0"\t"($4+$7)/$2*100"\t"($6/$3)*100}' | sort -k1,1 >> STAR_QC.tsv

    echo -e "Sample\tJunction_saturation_slope" > Junction_saturation.tsv
    tail -n +2 multiqc_data/multiqc_rseqc_junction_annotation.txt | cut -f1,21,20  | awk '{print $1"\t"($3-$2)/$3*100}' | sort -k1,1 >> Junction_saturation.tsv

    echo -e "Sample\tReads_mapping_sense_percentage\tReads_mapping_antisense_percentage\tReads_undetermined_strandedness_percentage" > strandedness_library_prep.tsv
    tail -n +2 multiqc_data/multiqc_rseqc_infer_experiment.txt | cut -f 1,2,3,4 | sort -k1,1 >> strandedness_library_prep.tsv

    awk '{print $1"\t"$2"\t"$3"\t"$4}' correlation_coefs.tsv | head -n 1 > new_corr.tsv
    awk '{print $1"\t"$2"\t"$3"\t"$4}' correlation_coefs.tsv | tail -n+2 | sort -k1,1 >> new_corr.tsv

    join -t $'\t' -j 1  samplenames.tsv  trackhub/trackhub_links.tsv | \\
    join -t $'\t' -j 1 -  fastqc_QC.tsv | \\
    join -t $'\t' -j 1 -  cutadapt_QC.tsv  | \\
    join -t $'\t' -j 1 -  STAR_QC.tsv  | \\
    join -t $'\t' -j 1 -  UMI_dedup_grouped.tsv  | \\
    join -t $'\t' -j 1 -  qualimap_QC.tsv  | \\
    join -t $'\t' -j 1 -  splicedReads_grouped.stats.tsv | \\
    join -t $'\t' -j 1 -  spliceJunctions_grouped.stats.tsv | \\
    join -t $'\t' -j 1 -  Junction_saturation.tsv | \\
    join -t $'\t' -j 1 -  read_coverage_uniformity_score.tsv | \\
    join -t $'\t' -j 1 -  insert_size.tsv | \\
    join -t $'\t' -j 1 -  library_balance.tsv | \\
    join -t $'\t' -j 1 -  strandedness_library_prep.tsv  | \\
    join -t $'\t' -j 1 -  biotype_table.tsv  | \\
    join -t $'\t' -j 1 -  new_corr.tsv > QC_table.tsv

    '''

}
