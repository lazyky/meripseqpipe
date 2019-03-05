#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : THREAD_NUM
Aligner_name=$1
THREAD_NUM=$2

for bam_file in *${Aligner_name}.bam
do
{
    samtools sort -@ ${THREAD_NUM:=1} $bam_file ${bam_file/.bam/_sort}
    samtools index ${bam_file/.bam/_sort.bam}
}&
done