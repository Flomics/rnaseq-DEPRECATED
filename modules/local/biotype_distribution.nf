process BIOTYPE_DISTRIBUTION {
    tag "$meta.id"
    label 'process_low'

    container "flomicsbiotech/flomics_qc_rnaseq:dev"

    input:
    tuple val(meta), path(bam)
    path  gtf

    output:
    tuple val(meta), path("*.tsv"), emit: results

    when:
    task.ext.when == null || task.ext.when

    shell:
    '''
    # read distribution by biotype

    # get exon records from gencode gtf
    awk '$3 == "exon"' !{gtf} > filtered_annotation_exon.gtf

    #get list of biotypes present in annotation
    cat filtered_annotation_exon.gtf | extractGffAttributeValue.pl gene_type|sort|uniq > biotypes_list.txt

    uuid=$(uuidgen)
    #make filtered BAM
    mkdir !{meta.id}
    cp !{bam} !{meta.id}
    samtools view -H !{meta.id}/!{bam} > /tmp/$uuid.sam
    samtools view -F2304 S!{meta.id}/!{bam} >> /tmp/$uuid.sam
    samtools view -b /tmp/$uuid.sam > !{meta.id}/!{meta.id}.umi_dedup.sorted.-F2304.bam

    #verify that we can select one and only one primary alignment per fragment
    if [ $(samtools view !{meta.id}/!{meta.id}.umi_dedup.sorted.-F2304.bam |cut -f1|sort|uniq -dc|awk '$1>2'|wc -l) -gt 0 ]; then
    echo "Found more than two primary alignments for some fragments. Cannot continue."
    exit 1
    fi

    bedtools intersect -split -abam !{meta.id}/!{meta.id}.umi_dedup.sorted.-F2304.bam -b filtered_annotation_exon.gtf -wo -bed | perl -F'\t' -slane 'chomp; $F[3]=~s/\/(?:1|2)$//; if($_=~/gene_type "(\S+)";/) {print "$F[3]\t$1"} else{die "No gene_type attribute found, cannot continue."}' | sort|uniq > tmp

    biotype_distribution.pl biotypes_list.txt tmp > biotypes_distribution.tsv
    '''
}
