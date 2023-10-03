#!/bin/bash
bam=$1
prefix=$2
project=$3
uuid=$4
profile=$5
s3_bucket_name="s3://flomics-public/RNAseq_pipeline/$project/$uuid/"
http_folder="https://flomics-public.s3.eu-west-1.amazonaws.com/RNAseq_pipeline/$project/$uuid/"

#Create the index files if they don't exist
for file in *.bam
        do if ! test -f ${file}.bai
        then samtools index $file
        fi
done

#Create hub.txt
mkdir -p trackhub
echo -e "hub $project\nshortLabel $project project\nlongLabel samples from the $project project\ngenomesFile genomes.txt\nemail lluc.cabus@flomics.com" > trackhub/hub.txt
echo -e "track $prefix\nbigDataUrl ${http_folder}dataFiles/$bam\nshortLabel $prefix\nlongLabel $prefix sample of the project $project\ntype bam\nvisibility full\ndoWiggle on\nmaxHeightPixels 100:32:10\n" > trackhub/${prefix}_trackDb.txt
echo "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hubUrl=${http_folder}hub.txt" > trackhub/UCSC.txt
echo -e "genome hg38\ntrackDb hg38/trackDb.txt" > trackhub/genomes.txt

#Create assembly hub
mkdir -p assembly_hub
echo -e "hub $project\nshortLabel $project project\nlongLabel samples from the $project project\ngenomesFile genomes.txt\nemail marc.weber@flomics.com" > assembly_hub/hub.txt
echo -e "track $prefix\nbigDataUrl ${http_folder}dataFiles/$bam\nshortLabel $prefix\nlongLabel $prefix sample of the project $project\ntype bam\nvisibility full\ndoWiggle on\nmaxHeightPixels 100:32:10\n" > assembly_hub/${prefix}_trackDb.txt
echo "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hubUrl=${http_folder}assembly_hub.hub.txt" >> assembly_hub/UCSC.txt
echo -e "genome hg38_ERCC_spikeins\ntrackDb hg38_ERCC_spikeins/trackDb.txt\ntwoBitPath hg38_ERCC_spikeins/hg38_ERCC_spikeins.2bit\norganism H. sapiens\ndefaultPos ERCC-00130:1-1051" > assembly_hub/genomes.txt

#Upload files to s3
if [[ $profile != *"test"* ]]; then
    aws s3 cp . ${s3_bucket_name}dataFiles/ --recursive --exclude "*" --include "*.bam*" --acl public-read
    aws s3 cp trackhub/hub.txt ${s3_bucket_name}trackhub/ --acl public-read
    aws s3 cp trackhub/genomes.txt ${s3_bucket_name}trackhub/ --acl public-read
    aws s3 cp assembly_hub/hub.txt ${s3_bucket_name}assembly_hub/ --acl public-read
    aws s3 cp assembly_hub/genomes.txt ${s3_bucket_name}assembly_hub/ --acl public-read
fi
