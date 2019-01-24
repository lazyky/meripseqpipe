#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : designfile
#$3 argv 3 : gtf file
#$4 argv 4 : THREAD_NUM
Aligner_name=$1
designfile=$2
gtf_file=$3
THREAD_NUM=$4
MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
lable=$(for ((i=0;i<=$MAX_SITUATION;i++));do echo "Situation"$i","; done |awk 'BEGIN{FS=",";OFS=",";ORS=" "}{print $1,$2}')
bam_file=$(for ((i=0;i<=$MAX_SITUATION;i++));do ls *input_${i}_${Aligner_name}*.bam | awk '{ORS=","}{print $0}END{ORS=" ";print "\t"}';done | awk 'BEGIN{FS=",\t";OFS=" "}{print $1,$2}')
echo $lable
echo $bam_file 
#cuffdiff -L $lable -p ${THREAD_NUM} --time-series --multi-read-correct --library-type fr-unstranded --poisson-dispersion $gtf_file ${bam_file}
