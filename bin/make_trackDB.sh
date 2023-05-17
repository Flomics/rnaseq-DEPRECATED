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
echo -e "hub $project\nshortLabel $project project\nlongLabel samples from the $project project\ngenomesFile genomes.txt\nemail lluc.cabus@flomics.com" > hub.txt
echo -e "track $prefix\nbigDataUrl ${http_folder}dataFiles/$bam\nshortLabel $prefix\nlongLabel $prefix sample of the project $project\ntype bam\nvisibility full\ndoWiggle on\nmaxHeightPixels 100:32:10\n" > ${prefix}_trackDb.txt
echo "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hubUrl=${http_folder}hub.txt" > UCSC.txt
echo -e "genome hg38\ntrackDb hg38/trackDb.txt" > genomes.txt

#Upload files to s3
if [[ $profile != *"test"* ]]; then
    aws s3 cp . ${s3_bucket_name}dataFiles/ --recursive --exclude "*" --include "*.bam*" --acl public-read
    aws s3 cp hub.txt $s3_bucket_name --acl public-read
    aws s3 cp genomes.txt $s3_bucket_name --acl public-read
fi
