#!/bin/bash
#bash htseq_count.sh <gtf> <strand_info> <THREAD_NUM>
#$1 argv 1 : gtf file
#$2 argv 2 : strand_info
#$3 argv 3 : THREAD_NUM
gtf_file=$1
strand_info=$2
THREAD_NUM=$3
## Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

for bam_file in *.input*.bam
do
read -u 9
{
    featureCounts -p -s $strand_info -a ${gtf_file} -o ${bam_file/%_sort*/}.txt  ${bam_file}
    echo >&9
}& 
done
wait
echo "done"
exec 9<&-
exec 9>&-
