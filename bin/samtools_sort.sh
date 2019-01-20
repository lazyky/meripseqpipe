#!/bin/bash
#$1 argv 1 : uesd Aligner
for bam_file in *$1.bam
do
    samtools sort $bam_file -o ${bam_file/.bam/_sort.bam}
    samtools index ${bam_file/.bam/_sort.bam}
done