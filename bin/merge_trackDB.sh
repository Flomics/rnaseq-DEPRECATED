#!/bin/bash

outdir=$1
timestamp=$2
project=$3
timestamp_simple=`echo "$timestamp" | sed 's/\..*$//' | sed 's/:/_/g'`

s3_bucket_name="s3://flomics-public/RNAseq_pipeline_trackHubs/$timestamp_simple/"

#Aggregate track_Db.txt files and uploads them to s3
for file in *trackDb.txt
do cat $file >> trackDb.txt
done
aws s3 cp trackDb.txt ${s3_bucket_name}hg38/ --acl public-read
