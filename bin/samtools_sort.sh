#!/bin/bash
#$1 argv 1 : uesd Aligner
Aligner_name=$1
for bam_file in *${Aligner_name}.bam
do
    samtools sort $bam_file ${bam_file/.bam/_sort}
    samtools index ${bam_file/.bam/_sort.bam}
done