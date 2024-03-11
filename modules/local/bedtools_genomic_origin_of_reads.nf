process BEDTOOLS_GENOMIC_ORIGIN_OF_READS {
    tag "$meta.id"
    label 'process_medium'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.txt"), emit: results
    tuple val(meta), path("*_genomic_origin_of_reads.tsv"), emit: table
    tuple val(meta), path("*_mqc.yaml"), emit: mqc_bedtools_goor_yaml
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    set -e

    #take only exon records
    awk '$3 == "exon"' !{gtf} > filtered_annotation_exon.gtf

    #extract gene locus coordinates
    cat !{gtf} | extract_locus_coords.pl - > gencode.genes.bed

    uuid=$(cat /proc/sys/kernel/random/uuid)
    #make filtered BAM
    mkdir !{meta.id}
    cp !{bam} !{meta.id}
    samtools view -H !{meta.id}/!{bam} > /tmp/$uuid.sam
    samtools view -F3332 !{meta.id}/!{bam} >> /tmp/$uuid.sam
    samtools view -b /tmp/$uuid.sam > !{meta.id}/!{meta.id}.umi_dedup.sorted.-F3332.bam

    #mapped fragments
    samtools view !{meta.id}/!{meta.id}.umi_dedup.sorted.-F3332.bam | cut -f1 | sort | uniq > !{meta.id}_mapped_fragments.list.txt

    #exonic fragments.
    #Do not take strandedness into account for simplicity. If we did, we would need to take the library type into account (e.g. only read#2 should match the strand in SMARTEr libraries)
    #Use -split to avoid taking into account gaps in spliced alignments
    # takes about 1.3Gb of mem without the -sorted option
    bedtools intersect -split -abam !{meta.id}/!{meta.id}.umi_dedup.sorted.-F3332.bam -b filtered_annotation_exon.gtf -wa -u  | samtools view - | cut -f1 | sort | uniq > !{meta.id}_exonic_fragments.list.txt

    #genic fragments
    bedtools intersect -split -abam !{meta.id}/!{meta.id}.umi_dedup.sorted.-F3332.bam -b gencode.genes.bed -wa -u  | samtools view - | cut -f1 | sort | uniq > !{meta.id}_genic_fragments.list.txt

    #intergenic fragments (mapped reads that do not overlap genic regions)
    comm -3 !{meta.id}_mapped_fragments.list.txt !{meta.id}_genic_fragments.list.txt > !{meta.id}_intergenic_fragments.list.txt

    #intronic fragments (genic fragments that do not overlap exons)
    comm -3 !{meta.id}_genic_fragments.list.txt !{meta.id}_exonic_fragments.list.txt  > !{meta.id}_intronic_fragments.list.txt

    #create the header for the csv output file
    echo "Exonic\tIntronic\tIntergenic\texonic_percentage\tintronic_percentage\tintergenic_percentage\tmapped_fragments" > !{meta.id}_genomic_origin_of_reads.tsv

    #variables for each field
    sample_name=!{meta.id}
    exonic_count=$(wc -l < !{meta.id}_exonic_fragments.list.txt)
    intronic_count=$(wc -l < !{meta.id}_intronic_fragments.list.txt)
    intergenic_count=$(wc -l < !{meta.id}_intergenic_fragments.list.txt)
    total=$((exonic_count+intronic_count+intergenic_count))
    exonic_percentage=$(echo "scale=2; ($exonic_count / $total) * 100" | bc)
    intronic_percentage=$(echo "scale=2; ($intronic_count / $total) * 100" | bc)
    intergenic_percentage=$(echo "scale=2; ($intergenic_count / $total) * 100" | bc)
    mapped_fragments=$(wc -l < !{meta.id}_mapped_fragments.list.txt)

    #populate csv
    echo "$exonic_count\t$intronic_count\t$intergenic_count\t$exonic_percentage\t$intronic_percentage\t$intergenic_percentage\t$mapped_fragments" >> !{meta.id}_genomic_origin_of_reads.tsv

    #create yaml for MultiQC
    echo -e "$sample_name: {Exonic: $exonic_count, Intronic: $intronic_count, Intergenic: $intergenic_count}" > !{meta.id}_genomic_origin_of_reads_mqc.yaml

    if [ $total == $mapped_fragments ]
    then
        echo "The sum of exonic, intronic, and intergenic fragments is equal to the number of mapped fragments."
    else
        echo "The sum of exonic, intronic, and intergenic fragments is NOT equal to the number of mapped fragments. Exiting."
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    !{task.process}:
        bedtools: $(bedtools --version | sed -e "s/bedtools v//g")
        samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
    END_VERSIONS
    '''
}
