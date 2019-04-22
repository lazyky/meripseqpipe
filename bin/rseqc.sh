#!/bin/bash
#$1 argv 1 : bed12
#$2 argv 2 : THREAD_NUM
bed12_file=$1
THREAD_NUM=$2

#Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

for bam_rseqc in *.bam
do 
read -u 9
{
    infer_experiment.py -i $bam_rseqc -r ${bed12_file} > ${bam_rseqc/_sort.bam/}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r ${bed12_file}
    bam_stat.py -i $bam_rseqc > ${bam_rseqc/_sort.bam/}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r ${bed12_file} 2> ${bam_rseqc/_sort.bam/}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r ${bed12_file}
    read_distribution.py -i $bam_rseqc -r ${bed12_file} > ${bam_rseqc/_sort.bam/}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.read_duplication
    echo >&9
}&
done
wait
echo "rseqc done"