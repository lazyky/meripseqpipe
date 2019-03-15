#!/bin/bash
#$1 argv 1 : fasta
#$2 argv 2 : gtf file
#$3 argv 3 : THREAD_NUM
fasta=$1
gtf_file=$2
THREAD_NUM=$3
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

for bed_file in *.bed
do
read -u 9
{
    annotatePeaks.pl ${bed_file} ${fasta} -gtf ${gtf_file} > ${bed_file/.bed/_annotatedbyhomer.bed}
    perl m6A_annotate_forGTF_xingyang2.pl ${gtf_file} ${bed_file} ${bed_file/.bed/_annotatedbyxy} 
    echo >&9
}& 
done
wait
echo "done"
exec 9<&-
exec 9>&-
