#!/bin/bash
#$1 argv 1 : fasta file
#$2 argv 2 : gtf file
#$3 argv 3 : THREAD_NUM
fasta_file=$1
gtf_file=$2
THREAD_NUM=$3
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

macs2_count=$(ls macs2_group*.summits| wc -w)
if [ $macs2_count -gt 1 ]; then
    for macs2_summits in macs2_group*.summits
    do
    read -u 9
    {
        sort -k5,5 -n -r ${macs2_summits} | head -1000 |awk '{summit=$3; print $1"\t"summit-51"\t"summit+50}' > ${macs2_summits/.summits/.location}
        intersectBed -wo -a ${macs2_summits/.summits/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${macs2_summits/.summits/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${macs2_summits/.summits/_bestpeaks.bed} -fo ${macs2_summits/.summits/_bestpeaks.fa}
        dreme -oc ${macs2_summits/.summits/_dreme} -p ${macs2_summits/.summits/_bestpeaks.fa} -rna
        echo >&9
    }& 
    done
fi

metpeak_count=$(ls metpeak_group*.bed| wc -w)
if [ $metpeak_count -gt 1 ]; then
    for metpeak_bed in metpeak_group*.bed
    do
    read -u 9
    {
        sort -k5,5 -n -r ${metpeak_bed} | head -1000 |awk '{ print $1"\t"$2"\t"$3}' > ${metpeak_bed/.bed/.location}
        intersectBed -wo -a ${metpeak_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${metpeak_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${metpeak_bed/.bed/_bestpeaks.bed} -fo ${metpeak_bed/.bed/_bestpeaks.fa}
        dreme -oc ${metpeak_bed/.bed/_dreme} -p ${metpeak_bed/.bed/_bestpeaks.fa} -rna
        echo >&9
    }& 
    done
fi

matk_count=$(ls MATK_group*.bed| wc -w)
if [ $matk_count -gt 1 ]; then
    for matk_bed in MATK_group*.bed
    do
    read -u 9
    {
        sort -k5,5 -n ${matk_bed} | head -1000 |awk '{ print $1"\t"$2"\t"$3}' > ${matk_bed/.bed/.location}
        intersectBed -wo -a ${matk_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${matk_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${matk_bed/.bed/_bestpeaks.bed} -fo ${matk_bed/.bed/_bestpeaks.fa}
        dreme -oc ${matk_bed/.bed/_dreme} -p ${matk_bed/.bed/_bestpeaks.fa} -rna
        echo >&9
    }& 
    done
fi

bedtools_count=$(ls bedtools_group*.bed| wc -w)
if [ $bedtools_count -gt 1 ]; then
    for bedtools_bed in bedtools_group*.bed
    do
    read -u 9
    {
        awk '{ print $1"\t"$2"\t"$3}' ${bedtools_bed} > ${bedtools_bed/.bed/.location}
        intersectBed -wo -a ${bedtools_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${bedtools_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${bedtools_bed/.bed/_bestpeaks.bed} -fo ${bedtools_bed/.bed/_bestpeaks.fa}
        dreme -oc ${bedtools_bed/.bed/_dreme} -p ${bedtools_bed/.bed/_bestpeaks.fa} -rna
        echo >&9
    }& 
    done
fi

mspc_count=$(ls mspc_group*.bed| wc -w)
if [ $mspc_count -gt 1 ]; then
    for mspc_bed in mspc_group*.bed
    do
    read -u 9
    {
        awk '{ print $1"\t"$2"\t"$3}' ${mspc_bed} > ${mspc_bed/.bed/.location}
        intersectBed -wo -a ${mspc_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${mspc_bed/.bed/_bestpeaks.bed}
        mspc getfasta -s -fi ${fasta_file} -bed ${mspc_bed/.bed/_bestpeaks.bed} -fo ${mspc_bed/.bed/_bestpeaks.fa}
        dreme -oc ${mspc_bed/.bed/_dreme} -p ${mspc_bed/.bed/_bestpeaks.fa} -rna
        echo >&9
    }& 
    done
fi

wait
echo "DREME done"
exec 9<&-
exec 9>&-