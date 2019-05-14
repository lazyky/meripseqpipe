#!/bin/bash
#bash cufflinks.sh <designfile> <gtf> <THREAD_NUM>
#$1 argv 1 : designfile
#$2 argv 2 : gtf file
#$3 argv 3 : THREAD_NUM
designfile=$1
gtf_file=$2
THREAD_NUM=$3

group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
tag=$(echo $group_list | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x"," }END{print $x" "}')

## Generate the array of bam files' name by groups eg. group1_1,group1_2 group2_1,group2_2 group3_1
bam_file_array=$(for group_id in $group_list 
        do
               echo *input_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x"," }END{print $x" "}'
        done 
        )

## Generate the information of assembly for cuffdiff 
rm -f assembly_list.txt
for bam_file in *input*.bam
do
        cufflinks -p ${THREAD_NUM} -G ${gtf_file} --library-type fr-unstranded -o ${bam_file/.bam/} $bam_file
        echo "./"${bam_file/.bam/}"/transcripts.gtf" >> assembly_list.txt
done
cuffmerge -o ./merged_gtf -g ${gtf_file} -p ${THREAD_NUM} assembly_list.txt

## Run Cuffdiff for differential expression analysis
cuffdiff -o cuffdiff\
         -L $tag \
         -p ${THREAD_NUM} \
         --time-series --multi-read-correct \
         --library-type fr-unstranded \
         ./merged_gtf/merged.gtf ${bam_file_array}
