#!/bin/bash

outdir=$1
project_name=${outdir#"s3://flomics-data/"}
folder_name=${project_name%/*}
folder_name=${folder_name%/*}
project_name=${project_name%%/*}

s3_bucket_name="s3://flomics-public/RNAseq_pipeline_trackHubs/"$folder_name"/"
data_file_folder=${s3_bucket_name}"dataFiles/"
#Aggregate track_Db.txt files and uploads them to s3
for file in *trackDb.txt
do cat $file >> trackDb.txt
done
aws s3 cp trackDb.txt ${s3_bucket_name}hg38/ --acl public-read
