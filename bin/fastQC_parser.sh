#!/bin/bash
num_ends_sequencing=$1
prefix=$2
cat multiqc_data/multiqc_fastqc_1.txt | grep "$prefix" | cut -f12,13,14,15,16,17,18,19,20 > tmp_fastQC.tsv
if [[ "$num_ends_sequencing" == "paired_end" ]]
then
  	for i in {1..9};
	do
          	row1=$(head -n1 tmp_fastQC.tsv| cut -f $i)
           	row2=$(tail -n1 tmp_fastQC.tsv| cut -f $i)
           	if [ "$row1" = "fail" ] || [ "$row2" = "fail" ]; then
                   	echo -ne "fail\t" >> ${prefix}_fastqc_QC.tsv
           	elif [ "$row1" == "warn" ] || [ "$row2" = "warn" ]; then
                   	echo -ne "warn\t" >> ${prefix}_fastqc_QC.tsv
           	else
                   	echo -ne "pass\t" >> ${prefix}_fastqc_QC.tsv
           	fi
    done
else
    cat tmp_fastQC.tsv >> ${prefix}_fastqc_QC.tsv
fi
sed -i 's/[ \t]*$//' ${prefix}_fastqc_QC.tsv
printf '\n' >> ${prefix}_fastqc_QC.tsv
