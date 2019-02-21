#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : gtf file
#$3 argv 3 : THREAD_NUM
Aligner_name=$1
gtf_file=$2
THREAD_NUM=$3
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

for bam_file in *${Aligner_name}*.bam
do
read -u 9
{
    htseq-count -f bam --stranded=no ${bam_file} ${gtf_file} > ${bam_file/%_${Aligner_name}_*/_${Aligner_name}}.txt
    echo >&9
}& 
done
wait
echo "done"
exec 9<&-
exec 9>&-
