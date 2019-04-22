#!/bin/bash
## $1 argv 1 : matk_jar
## $2 argv 2 : fasta file
## $3 argv 3 : gtf file
## $4 argv 4 : THREAD_NUM
matk_jar=$1
fasta_file=$2
gtf_file=$3
THREAD_NUM=$4

mkfifo tmp
exec 9<>tmp
#rm -rf /tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

#check if the output file of Macs2 exists
bed_count=$(ls *.bed |wc -w)
if [ $bed_count -gt 0 ]; then
    for bedfile in $(ls *.bed)
    do
    read -u 9
    {
        awk '{print $1"\t"$2"\t"$3}' ${bedfile} > ${bedfile/.bed/.location}
        intersectBed -wo -a ${bedfile/.bed/.location} -b ${gtf_file} | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${bedfile/.bed/_bestpeaks.bed}
        bedtools getfasta -s -fi ${fasta_file} -bed ${bedfile/.bed/_bestpeaks.bed} -fo ${bedfile/.bed/_bestpeaks.fa}
        java -jar $matk_jar -singleNucleotide -mode "Fasta" -fasta  ${bedfile/.bed/_bestpeaks.fa} -out ${bedfile/.bed/_matk_prediction}.txt
        perl runsramp.pl ${bedfile/.bed/_bestpeaks.fa} tmp.${bedfile/.bed/_sramp.txt} full
        head -1 tmp.${bedfile/.bed/_sramp.txt} > ${bedfile/.bed/_sramp_prediction.txt}
        cat tmp.${bedfile/.bed/_sramp.txt} | grep "m6A site (" >> ${bedfile/.bed/_sramp_prediction.txt}
        rm tmp*
        echo >&9
    }& 
    done
fi
wait
echo "m6Aprediction done"
exec 9<&-
exec 9>&-