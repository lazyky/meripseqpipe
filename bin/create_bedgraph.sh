#!/bin/bash
#$1 argv 1 : THREAD_NUM
THREAD_NUM=$1
for bam_file in *.bam
do
{
    bamCoverage -b $bam_file --outFileFormat bedgraph -o ${bam_file%.bam*}.bedgraph -p ${THREAD_NUM:=1}
}
done
