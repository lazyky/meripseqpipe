#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : gtf file
#$2 argv 3 : Read classification
Aligner_name=$1
gtf_file=$2
input_or_ip=$3
for bam_file in *$3*$1*.bam
do 
    htseq-count -f bam -s no $bam_file $2 > ${bam_file%_$1_*}_$1.txt
done
