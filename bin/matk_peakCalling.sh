#!/bin/bash
## MATK_peakCalling.sh tophat2 $matk_jar $designfile
## $1 argv 1 : uesd Aligner
## $2 argv 2 : matk_jar
## $3 argv 3 : designfile
Aligner_name=$1
matk_jar=$2
designfile=$3
MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
for ((i=0;i<=$MAX_SITUATION;i++))
do 
    ip_bam_file_array=$(echo *ip_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    input_bam_file_array=$(echo *input_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    java -jar $matk_jar -peakCalling -ip "$ip_bam_file_array" -input "$input_bam_file_array" -out MATK_peakCalling_situation_${i}.bed
done 
