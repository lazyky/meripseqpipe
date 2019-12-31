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
    chrom_size=$4
    prefix=$5
    sort -k5,5 -g ${bed_file}| head -1000 | awk '{ print $1"\t"$2"\t"$3}' > ${prefix}.location
    intersectBed -wo -a ${prefix}.location -b $gtf | awk -v OFS="\t" '{print $1,$2,$3,"*","*",$10}' | sort -k1,2 | uniq > ${prefix}_bestpeaks.bed
    ame -oc ${prefix}_ame ${prefix}_bestpeaks.fa $RRACH_motif
    fastaFromBed -name+ -split -s -fi $fasta -bed ${prefix}_bestpeaks.bed > ${prefix}_bestpeaks.fa
    shuffleBed -incl ${bed_file}| -seed 12345 -noOverlapping -i ${prefix}_bestpeaks.bed -g $chrom_size > ${prefix}_random_peak.bed
    fastaFromBed -name+ -split -s -fi $fasta -bed ${prefix}_random_peak.bed > ${prefix}_random_peak.fa
    findMotifs.pl ${prefix}_bestpeaks.fa fasta ${prefix}_homer -fasta ${prefix}_random_peak.fa -p ${THREAD_NUM:=1} \
        -len 5,6,7,8 -S 10 -rna -dumpFasta > ${prefix}_homer_run.log 2>&1
    #dreme -k 7 -oc ${prefix}_dreme -p ${prefix}_bestpeaks.fa -rna
}

#check if the output file of Bedtools Merge exists
bed_count=$(ls *.bed| wc -w)
cat $fasta_file | awk 'BEGIN{len=""}{if($0~">"){split($0,ID,"[> ]");printf len"ABC"ID[2]"\t";len=0}else{len=len+length($0)}}END{print len}' |sed 's/ABC/\n/g' |awk NF > chromsizes.file
if [ $bed_count -gt 0 ]; then
    for bedfile in *.bed
    do
    {
        motif_searching_by_pvalue $bedfile $fasta_file $gtf_file chromsizes.file ${bedfile/.bed/}
    }
    done
fi
wait
echo "Searching motif done"
exec 9<&-
exec 9>&-