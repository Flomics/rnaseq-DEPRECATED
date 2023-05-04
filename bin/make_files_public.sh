#!/bin/bash
multiqc_report=$1
flomics_QC_report_table=$2
project=$3
uuid=$4
s3_bucket_name="s3://flomics-public/RNAseq_pipeline/$project/$uuid/"

aws s3 cp $multiqc_report $s3_bucket_name --acl public-read
aws s3 cp $flomics_QC_report_table $s3_bucket_name --acl public-read