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

for macs2_summits in macs2_situation*.summits
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

for metpeak_bed in metpeak_situation*.bed
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

wait
echo "DREME done"
exec 9<&-
exec 9>&-
java -jar /data2/yeying_by_zky/EBV_pipe/bin/MATK-1.0.jar -diff
                         -control_ip RE-C2E-IP.read1_Clean_control_ip_1_star_sort.bam                         
                         -control_input RE-C2E-input.read1_Clean_control_input_1_star_sort.bam                        
                          -treated_ip RE-C2EN-IP.read1_Clean_treated_ip_2_star_sort.bam                        
                           -treated_input RE-C2EN-input.read1_Clean_treated_input_2_star_sort.bam                         
                           -control_bed MATK_situation_1_star.bed                         
                           -treated_bed MATK_situation_2_star.bed                         
                           -gtf EBV_Akata.gtf                         
                           -out diffMATK_situation_1__2_star.bed