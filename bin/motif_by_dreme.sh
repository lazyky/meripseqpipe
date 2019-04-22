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

#check if the output file of Macs2 exists
bed_count=$(ls *.bed| grep -v bedtools |wc -w)
if [ $bed_count -gt 0 ]; then
    for bedfile in $(ls *.bed| grep -v bedtools| grep -v mspc)
    do
    read -u 9
    {
        sort -k5,5 -n -r ${bedfile} | head -1000 |awk '{ print $1"\t"$2"\t"$3}' > ${bedfile/.bed/.location}
        intersectBed -wo -a ${bedfile/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${bedfile/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${bedfile/.bed/_bestpeaks.bed} -fo ${bedfile/.bed/_bestpeaks.fa}
        dreme -k 6 -oc ${bedfile/.bed/_dreme} -p ${bedfile/.bed/_bestpeaks.fa} -rna
        findMotifsGenome.pl ${bedfile/.bed/_bestpeaks.bed} ${fasta_file} ${bedfile/.bed/_homer} -len 6 -rna
        echo >&9
    }& 
    done
fi

#check if the output file of Bedtools Merge exists
bedtools_count=$(ls bedtools_group*.bed| wc -w)
if [ $bedtools_count -gt 0 ]; then
    for bedtools_bed in bedtools_group*.bed
    do
    read -u 9
    {
        sort -k5,5 -n -r ${bedtools_bed}| head -1000 | awk '{ print $1"\t"$2"\t"$3}' > ${bedtools_bed/.bed/.location}
        intersectBed -wo -a ${bedtools_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${bedtools_bed/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${bedtools_bed/.bed/_bestpeaks.bed} -fo ${bedtools_bed/.bed/_bestpeaks.fa}
        dreme -k 6 -oc ${bedtools_bed/.bed/_dreme} -p ${bedtools_bed/.bed/_bestpeaks.fa} -rna
        findMotifsGenome.pl ${bedtools_bed/.bed/_bestpeaks.bed} ${fasta_file} ${bedtools_bed/.bed/_homer} -len 6 -rna
        echo >&9
    }& 
    done
fi

#check if the output file of MSPC Merge exists
mspc_count=$(ls mspc_group*.bed| wc -w)
if [ $mspc_count -gt 0 ]; then
    for mspc_bed in mspc_group*.bed
    do
    read -u 9
    {
        sort -k5,5 -n -r ${mspc_bed} | head -1000 | awk '{ print $1"\t"$2"\t"$3}' > ${mspc_bed/.bed/.location}
        intersectBed -wo -a ${mspc_bed/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${mspc_bed/.bed/_bestpeaks.bed}
        mspc getfasta -s -fi ${fasta_file} -bed ${mspc_bed/.bed/_bestpeaks.bed} -fo ${mspc_bed/.bed/_bestpeaks.fa}
        dreme -k 6 -oc ${mspc_bed/.bed/_dreme} -p ${mspc_bed/.bed/_bestpeaks.fa} -rna
        findMotifsGenome.pl ${mspc_bed/.bed/_bestpeaks.bed} ${fasta_file} ${mspc_bed/.bed/_homer} -len 6 -rna
        echo >&9
    }& 
    done
fi

wait
echo "DREME done"
exec 9<&-
exec 9>&-