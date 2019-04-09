#!/bin/bash
#$1 argv 1 : gtf file
#$2 argv 2 : THREAD_NUM
gtf_file=$1
THREAD_NUM=$2
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

for bam_file in *.bam
do
read -u 9
{
    htseq-count -f bam --stranded=no ${bam_file} ${gtf_file} > ${bam_file/%_sort*/}.txt
    echo >&9
}& 
done
wait
echo "done"
exec 9<&-
exec 9>&-
