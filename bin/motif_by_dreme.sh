#!/bin/bash
#$1 argv 1 : fasta file
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

for macs2_bed in macs2_situation*.bed
do
read -u 9
{
    do
        sort -k5,5 -n -r ${macs2_bed} | head -1000 |awk '{summit=$3; print $1"\t"summit-51"\t"summit+50}' > ${macs2_bed/.bed/.location}
        intersectBed -wo -a ${macs2_bed/.bed/.location} -b !{gtf} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${macs2_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi !{fasta} -bed ${macs2_bed/.bed/_bestpeaks.bed} -fo ${macs2_bed/.bed/_bestpeaks.fa}
        dreme -oc ${macs2_bed/.bed/_dreme} -p ${macs2_bed/.bed/_bestpeaks.fa} -rna
    done
}& 
done

for metpeak_bed in metpeak_situation*.bed
do
read -u 9
{
    do
        sort -k5,5 -n -r ${metpeak_bed} | head -1000 |awk '{ print $1"\t"$2"\t"$3}' > ${metpeak_bed/.bed/.location}
        intersectBed -wo -a ${metpeak_bed/.bed/.location} -b !{gtf} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${macs2_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi !{fasta} -bed ${metpeak_bed/.bed/_bestpeaks.bed} -fo ${metpeak_bed/.bed/_bestpeaks.fa}
        dreme -oc ${metpeak_bed/.bed/_dreme} -p ${metpeak_bed/.bed/_bestpeaks.fa} -rna
    done
}& 
done

wait
echo "done"
exec 9<&-
exec 9>&-
