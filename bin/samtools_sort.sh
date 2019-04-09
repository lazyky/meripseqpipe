#!/bin/bash
#$1 argv 1 : THREAD_NUM
THREAD_NUM=$1

for bam_file in *.bam
do
{
    samtools sort -@ ${THREAD_NUM:=1} $bam_file ${bam_file/.bam/_sort}
    samtools index ${bam_file/.bam/_sort.bam}
}
done
wait
echo " Bam files are sorted"
