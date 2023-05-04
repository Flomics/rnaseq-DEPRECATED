process FLOMICS_QC_AGGREGATOR{
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:latest"


    input:
    path multiqc_data
    path trackhub_links
    path trackDbs
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


    '''
    merge_trackDB.sh !{project} !{uuid} #Merges the trackDb files and uploads it to s3

    echo -e "trackhub_link" > trackhub_links.tsv
    for file in *_trackhub_links.tsv
        do cat $file >> trackhub_links.tsv
    done

    biotype_table_parser.r  #Sorts the columns of the biotypes according to last table
    gene_coverage_profile_calculation.r #Calculates the gene coverage profile


    echo -e "per_base_sequence_quality\tper_sequence_quality_scores\tper_base_sequence_content\tper_sequence_gc_content\tper_base_n_content\tsequence_length_distribution\tsequence_duplication_levels\toverrepresented_sequences\tadapter_content" > fastqc_QC.tsv
    for file in *_fastqc_QC.tsv; do
        cat $file >> fastqc_QC.tsv
    done


    echo -e "uniquely_mapped_reads\tspliced_reads\tpercentage_of_spliced_reads" > splicedReads_grouped.stats.tsv
    for file in *splicedReads.stats.tsv; do
        tail -n +2 $file >> splicedReads_grouped.stats.tsv
    done

    echo -e "total_splice_junctions\tknown_splice_junctions\t%known_splice_junctions" > spliceJunctions_grouped.stats.tsv
    for file in *spliceJunctions.stats.tsv; do
        tail -n +2 $file >> spliceJunctions_grouped.stats.tsv
    done

    echo -e "median_insert_size" > insert_size.tsv
    for file in *.insert_size_median.tsv; do
        tail -n +2 $file >> insert_size.tsv
    done

    echo -e "number_of_genes_contributing_to_1%_of_reads\tnumber_of_genes_contributing_to_5%_of_reads\tnumber_of_genes_contributing_to_10%_of_reads\tnumber_of_genes_contributing_to_50%_of_reads\tnumber_of_genes_contributing_to_80%_of_reads" > library_balance.tsv

    for file in *_genes_contributing_to_percentage_reads.tsv; do
        tail -n +2 $file >> library_balance.tsv
    done

    echo -e "Number_mapped_unique_molecules\tPercentage_unique_molecules" > UMI_dedup_grouped.tsv
    if ls *_UMI_dedup.tsv 1> /dev/null 2>&1; then
        for file in *_UMI_dedup.tsv; do
            tail -n +2 $file >> UMI_dedup_grouped.tsv
        done
    else
        for file in *splicedReads.stats.tsv; do
            echo -e "NA\tNA" >> UMI_dedup_grouped.tsv
        done
    fi

    cut -f1 multiqc_data/mqc_qualimap_genomic_origin_1.txt > samplenames.tsv
    echo -e "Read_number\tReads_passing_trimming\tPercentage_reads_passing_trimming"> cutadapt_QC.tsv
    cut -f1 multiqc_data/multiqc_cutadapt.txt | sed 's/_[0-9]$//g' > tmp.cutadapt.txt
    paste tmp.cutadapt.txt multiqc_data/multiqc_cutadapt.txt | awk '!seen[$1]++' > tmp2.cutadapt.txt
    cut -f 4,6 tmp2.cutadapt.txt | tail -n +2 | awk '{print $1"\t"$2"\t"($2/$1*100)}' >> cutadapt_QC.tsv

    echo -e "Exonic\tIntronic\tIntergenic\tExonic_percentage\tIntronic_percentage\tIntergenic_percentage" > qualimap_QC.tsv
    tail -n +2 multiqc_data/mqc_qualimap_genomic_origin_1.txt | cut --complement -f1 | awk '{print $0"\t"$1/($1+$2+$2)*100"\t"$2/($1+$2+$2)*100"\t"$3/($1+$2+$2)*100'} >> qualimap_QC.tsv

    echo -e "total_reads\tavg_input_read_length\tnumber_of_uniquely_mapped_reads\tpercentage_of_uniquely_mapped_reads\tavg_mapped_read_length\tnumber_of_multimapped_reads\tpercentage_of_unmapped_too_short_reads\tmapped_percentage\taverage_mapped_length_percentage" > STAR_QC.tsv
    tail -n +2  multiqc_data/multiqc_star.txt | cut -f2,3,4,5,6,18,23 | awk '{print $0"\t"($3+$6)/$1*100"\t"($5/$2)*100}' >> STAR_QC.tsv

    echo "Junction_saturation_slope" > Junction_saturation.tsv
    tail -n +2 multiqc_data/mqc_rseqc_junction_saturation_plot_All_Junctions.txt | cut -f21,20  | awk '{print ($2-$1)/$2*100}' >> Junction_saturation.tsv

    echo -e "Reads_mapping_sense_percentage\tReads_mapping_antisense_percentage\tReads_undetermined_strandedness_percentage" > strandedness_library_prep.tsv
    tail -n +2 multiqc_data/mqc_rseqc_infer_experiment_plot_1.txt | cut -f 2,3,4 >> strandedness_library_prep.tsv

    awk '{print $2"\t"$3"\t"$4}' correlation_coefs.tsv > new_corr.tsv

    paste samplenames.tsv trackhub_links.tsv fastqc_QC.tsv cutadapt_QC.tsv STAR_QC.tsv UMI_dedup_grouped.tsv qualimap_QC.tsv splicedReads_grouped.stats.tsv \\
    spliceJunctions_grouped.stats.tsv Junction_saturation.tsv read_coverage_uniformity_score.tsv insert_size.tsv library_balance.tsv \\
    strandedness_library_prep.tsv biotype_table.tsv new_corr.tsv > QC_table.tsv

    '''



}
