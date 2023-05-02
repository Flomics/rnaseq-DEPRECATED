#!/bin/bash
multiqc_report=$1
flomics_QC_report_table=$2
timestamp=$3
project=$4
timestamp_simple=`echo "$timestamp" | sed 's/\..*$//' | sed 's/:/_/g'`
s3_bucket_name="s3://flomics-public/RNAseq_pipeline_trackHubs/${project}_$timestamp_simple/"

aws s3 cp $multiqc_report $s3_bucket_name --acl public-read
aws s3 cp $flomics_QC_report_table $s3_bucket_name --acl public-read
