#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : gtf file
for bam_file in *input*$1*.bam
do 
    htseq-count -f bam -s no $bam_file $2 > ${bam_file%_$1_*}_$1.txt
done
