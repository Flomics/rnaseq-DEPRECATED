#!/bin/bash
project=$1
uuid=$2
profile=$3
s3_bucket_name="s3://flomics-public/RNAseq_pipeline/$project/$uuid/"

#Aggregate track_Db.txt files and uploads them to s3
for file in trackhub/*_trackDB.txt
do
    cat $file >> trackhub/trackDb.txt
done

if [[ $profile != *"test"* ]]; then
    aws s3 cp trackhub/trackDb.txt ${s3_bucket_name}trackhub/hg38/ --acl public-read
fi

# Aggregate assembly hub track_Db.txt files and uploads them to s3
# We use a parent composite track (group of tracks) to simplify
# the configuration of all the tracks.
# See https://genome.ucsc.edu/goldenpath/help/hubQuickStartGroups.html#composite
cat << EOF > assembly_hub/trackDb.txt
track all_samples
compositeTrack on
shortLabel ${project}_all_samples
longLabel ${project}_all_samples
type bam
allButtonPair on
visibility full
autoscale on
maxHeightPixels 500:32:10
visibility full
doWiggle on


EOF
for file in assembly_hub/*_trackDb.txt
do
    cat $file >> assembly_hub/trackDb.txt
done

if [[ $profile != *"test"* ]]; then
    aws s3 cp assembly_hub/trackDb.txt ${s3_bucket_name}assembly_hub/hg38_ERCC_spikeins/ --acl public-read
    # We also copy the multi-region bed file from the genome reference folder to the
    # assembly hub folder on the public bucket
    aws s3 cp s3://flomics-no-backup/references/Genomes/Homo_sapiens/hg38_ERCC_spikeins/multiregion_all_spikeins.bed ${s3_bucket_name}assembly_hub/hg38_ERCC_spikeins/ --acl public-read
fi
