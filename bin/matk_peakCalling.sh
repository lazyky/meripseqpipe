#!/bin/bash
## MATK_peakCalling.sh $matk_jar $designfile
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : THREAD_NUM
## $4 argv 4 : flag_peakCallingbygroup
matk_jar=$1
designfile=$2
THREAD_NUM=$3
flag_peakCallingbygroup=$4

if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do 
    {
        ip_bam_file_array=$(echo *ip_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        input_bam_file_array=$(echo *input_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        java -jar $matk_jar -peakCalling -ip "$ip_bam_file_array" -input "$input_bam_file_array" -out MATK_group_${group_id}.bed
    }
    done 
else
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    {
        ip_bam_file=$(ls ${sample_id}.ip*.bam)
        input_bam_file=$(ls ${sample_id}.ip*.bam)
        java -jar $matk_jar -peakCalling -ip "$ip_bam_file" -input "$input_bam_file" -out MATK_${sample_id}.bed
    }
    done
fi
wait
echo "MATK done"