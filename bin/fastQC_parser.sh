#!/bin/bash
num_ends_sequencing=$1
END=$(wc -l <multiqc_data/multiqc_fastqc_1.txt)
if [[ "$num_ends_sequencing" == "paired_end" ]]
then
	cut -f11,12,13,14,15,16,17,18,19,20 multiqc_data/multiqc_fastqc_1.txt | head -n 1 > fastqc_QC.tsv
	CURRENT=2
	while [[ $CURRENT -le $END ]]; do
        	for i in {1..10}; do
                	row1=$(cut -f11,12,13,14,15,16,17,18,19,20 multiqc_data/multiqc_fastqc_1.txt | tail -n +$CURRENT | head -n 1 | cut -f $i)
                	row2=$(cut -f11,12,13,14,15,16,17,18,19,20 multiqc_data/multiqc_fastqc_1.txt | tail -n +$CURRENT | head -n 2 | tail -n 1 | cut -f $i)
                	if [ "$row1" = "fail" ] || [ "$row2" = "fail" ]; then
                        	echo -ne "fail\t" >> fastqc_QC.tsv
                	elif [ "$row1" == "warn" ] || [ "$row2" = "warn" ]; then
                        	echo -ne "warn\t" >> fastqc_QC.tsv
                	else
                        	echo -ne "pass\t" >> fastqc_QC.tsv
                	fi
        	done
        	echo >> fastqc_QC.tsv
        	CURRENT=$((CURRENT+2))
	done
else
        cut -f11,12,13,14,15,16,17,18,19,20 multiqc_data/multiqc_fastqc_1.txt > fastqc_QC.tsv 
fi
sed -i 's/[ \t]*$//' fastqc_QC.tsv
