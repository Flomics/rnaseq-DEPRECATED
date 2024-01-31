process BEDTOOLS_GENOMIC_ORIGIN_OF_READS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0' :
        'quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0' }"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.txt"), emit: results
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    #take only exon records
    awk '$3 == "exon"' !{gtf} > filtered_annotation_exon.gtf

    #name sort input bam
    samtools sort -n !{bam} -o name_sorted.bam

    #extract gene locus coordinates
    cat !{gtf} | extract_locus_coords.pl - > gencode.genes.bed

    #mapped fragments
    samtools view -F4 name_sorted.bam | cut -f1 | sort | uniq > "!{meta.id}_mapped_fragments.list.txt"

    #exonic fragments.
    #Do not take strandedness into account for simplicity. If we did, we would need to take the library type into account (e.g. only read#2 should match the strand in SMARTEr libraries)
    #Use -split to avoid taking into account gaps in spliced alignments
    # takes about 1.3Gb of mem without the -sorted option
    bedtools intersect -split -abam name_sorted.bam -b filtered_annotation_exon.gtf -wa -u  | samtools view -F4 - | cut -f1 | sort | uniq > "!{meta.id}_exonic_fragments.list.txt"

    #genic fragments
    bedtools intersect -split -abam name_sorted.bam -b gencode.genes.bed -wa -u  | samtools view -F4 - | cut -f1 | sort | uniq > "!{meta.id}_genic_fragments.list.txt"

    #intergenic fragments (mapped reads that do not overlap genic regions)
    comm -3 mapped_fragments.list.txt genic_fragments.list.txt > "!{meta.id}_intergenic_fragments.list.txt"

    #intronic fragments (genic fragments that do not overlap exons)
    comm -3 genic_fragments.list.txt exonic_fragments.list.txt  > "!{meta.id}_intronic_fragments.list.txt"

    '''

    // """
    // bedtools_genomic_origin_of_reads.sh $bam $gtf $meta.id

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    // END_VERSIONS
    // """
}
