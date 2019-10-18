#!/bin/bash
#bash motif_searching.sh <fasta> <gtf> <THREAD_NUM>
#$1 argv 1 : fasta file
#$2 argv 2 : gtf file
#$3 argv 3 : RRACH_motif file
#$4 argv 4 : THREAD_NUM
fasta_file=$1
gtf_file=$2
RRACH_motif=$3
THREAD_NUM=$4


## setting function for motif searching 
function motif_searching_by_pvalue()
{
    bed_file=$1
    fasta=$2
    gtf=$3
    prefix=$4
    sort -k5,5 -g ${bed_file}| head -1000 | awk '{ print $1"\t"$2"\t"$3}' > ${prefix}.location
    intersectBed -wo -a ${prefix}.location -b $gtf | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${prefix}_bestpeaks.bed
    bedtools getfasta -s -fi $fasta -bed ${prefix}_bestpeaks.bed -fo ${prefix}_bestpeaks.fa
    ame -oc ${prefix}_ame ${prefix}_bestpeaks.fa $RRACH_motif
    findMotifsGenome.pl ${prefix}_bestpeaks.bed $fasta ${prefix}_homer -len 7 -rna -p ${THREAD_NUM:=1}
    #dreme -k 7 -oc ${prefix}_dreme -p ${prefix}_bestpeaks.fa -rna
}

#check if the output file of Bedtools Merge exists
bed_count=$(ls *.bed| wc -w)
if [ $bed_count -gt 0 ]; then
    for bedfile in *.bed
    do
    {
        motif_searching_by_pvalue $bedfile $fasta_file $gtf_file ${bedfile/.bed/}
    }
    done
fi
wait
echo "Searching motif done"
exec 9<&-
exec 9>&-