#!/bin/bash
bam=$1
outdir=$2
prefix=$3
folder_name=${outdir#"s3://flomics-data/"}
project_name=${folder_name%%/*}

s3_bucket_name="s3://flomics-public/RNAseq_pipeline_trackHubs/test/"$folder_name
data_file_folder=${s3_bucket_name}"dataFiles/"
for file in *.bam
        do if ! test -f ${file}.bai
        then samtools index $file
        fi
done

aws s3 cp . $data_file_folder --recursive --exclude "*" --include "*.bam*" --acl public-read #Upload files and make them public

#Create hub.txt
echo -e "hub $project_name\nshortLabel $project_name project\nlongLabel samples from the $project_name project\ngenomesFile genomes.txt\nemail lluc.cabus@flomics.com" > hub.txt
echo -e "track $prefix\nbigDataUrl https://flomics-public.s3.eu-west-1.amazonaws.com/RNAseq_pipeline_trackHubs/${folder_name}dataFiles/$bam\nshortLabel $sample_name\nlongLabel $sample_name sample of the project $project_name\ntype bam\nvisibility full\ndoWiggle on\n" > ${prefix}_trackDb.txt
echo "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?hubUrl=https://flomics-public.s3.eu-west-1.amazonaws.com/RNAseq_pipeline_trackHubs/test/${folder_name}hub.txt" > UCSC.txt
echo -e "genome hg38\ntrackDb hg38/trackDb.txt" > genomes.txt

aws s3 cp hub.txt $s3_bucket_name --acl public-read
aws s3 cp genomes.txt $s3_bucket_name --acl public-read
