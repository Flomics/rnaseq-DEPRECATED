#!/bin/bash
project=$1
uuid=$2
s3_bucket_name="s3://flomics-public/RNAseq_pipeline/$project/$uuid/"

#Aggregate track_Db.txt files and uploads them to s3
for file in *trackDb.txt
do cat $file >> trackDb.txt
done
aws s3 cp trackDb.txt ${s3_bucket_name}hg38/ --acl public-read
