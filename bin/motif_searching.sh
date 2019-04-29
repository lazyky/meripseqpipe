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
function motif_searching_by_pvalue()
{
    bed_file=$1
    fasta=$2
    gtf=$3
    sort -k5,5 -n -r ${bed_file}| head -1000 | awk '{ print $1"\t"$2"\t"$3}' > ${bed_file/.bed/.location}
    intersectBed -wo -a ${bed_file/.bed/.location} -b $gtf | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${bed_file/.bed/_bestpeaks.bed}
    bedtools getfasta -s -fi $fasta -bed ${bed_file/.bed/_bestpeaks.bed} -fo ${bed_file/.bed/_bestpeaks.fa}
    dreme -k 6 -oc ${bed_file/.bed/_dreme} -p ${bed_file/.bed/_bestpeaks.fa} -rna
    findMotifsGenome.pl ${bed_file/.bed/_bestpeaks.bed} $fasta ${bed_file/.bed/_homer} -len 6 -rna
}

#check if the output file of peakCalling processes exists
# bed_count=$(ls *.bed| grep -v bedtools |wc -w)
# if [ $bed_count -gt 0 ]; then
#     for bedfile in $(ls *.bed| grep -v bedtools| grep -v mspc)
#     do
#     read -u 9
#     {
#         motif_searching_by_pvalue $bedfile $fasta_file $gtf_file
#         echo >&9
#     }& 
#     done
# fi

#check if the output file of Bedtools Merge exists
bedtools_count=$(ls bedtools_group*.bed| wc -w)
if [ $bedtools_count -gt 0 ]; then
    for bedtools_bed in bedtools_group*.bed
    do
    read -u 9
    {
        motif_searching_by_pvalue $bedfile $fasta_file $gtf_file
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
        motif_searching_by_pvalue $bedfile $fasta_file $gtf_file
        echo >&9
    }& 
    done
fi

wait
echo "DREME done"
exec 9<&-
exec 9>&-