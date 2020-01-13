#!/bin/bash
#bash geneBody_coverage2.sh <gtf> <THREAD_NUM>
#$1 argv 1 : gtf file
#$2 argv 2 : THREAD_NUM
bed12_file=$1
THREAD_NUM=$2
## Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9   
done

for bigwig_file in *.bigwig
do
read -u 9
{
    geneBody_coverage2.py -i $bigwig_file -o ${bigwig_file%.bigwig*}.rseqc.txt -r ${bed12_file}
    echo >&9
}& 
done
wait
echo "Calculate coverage of data is finish"