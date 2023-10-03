#!/bin/bash
project=$1
uuid=$2
profile=$3
s3_bucket_name="s3://flomics-public/RNAseq_pipeline/$project/$uuid/"

#Aggregate track_Db.txt files and uploads them to s3
for file in trackhub/*_trackDB.txt
do cat $file >> trackhub/trackDb.txt
done

if [[ $profile != *"test"* ]]; then
    aws s3 cp trackhub/trackDb.txt ${s3_bucket_name}trackhub/hg38/ --acl public-read
fi

#Aggregate assembly hub track_Db.txt files and uploads them to s3
for file in assembly_hub/*_trackDb.txt
do cat $file >> assembly_hub/trackDb.txt
done

if [[ $profile != *"test"* ]]; then
    aws s3 cp assembly_hub/trackDb.txt ${s3_bucket_name}assembly_hub/hg38_ERCC_spikeins/ --acl public-read
fi
