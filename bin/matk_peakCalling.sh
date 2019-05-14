#!/bin/bash
## MATK_peakCalling.sh <matk_jar> <designfile> <flag_peakCallingbygroup>
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : flag_peakCallingbygroup
## flag_peakCallingbygroup: 1(group) 0(sample)
matk_jar=$1
designfile=$2
flag_peakCallingbygroup=$3

# Check the mode of peakcalling
if [ $flag_peakCallingbygroup -gt 0 ]; then
    ## Running the peakcalling mode of MATK per group
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do 
    {
        ip_bam_file_array=$(echo *ip_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        input_bam_file_array=$(echo *input_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        java -jar $matk_jar -peakCalling -ip "$ip_bam_file_array" -input "$input_bam_file_array" -out MATK_group_${group_id}.bed
        awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,-log($5)/log(10)}' MATK_group_${group_id}.bed > MATK_group_${group_id}_normalized.bed
    }
    done 
else
    ## Running the peakcalling mode of MATK per sample
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    {
        ip_bam_file=$(ls ${sample_id}.ip*.bam)
        input_bam_file=$(ls ${sample_id}.input*.bam)
        java -jar $matk_jar -peakCalling -ip "$ip_bam_file" -input "$input_bam_file" -out MATK_${sample_id}.bed
        awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,-log($5)/log(10)}' MATK_${sample_id}.bed > MATK_${sample_id}_normalized.bed
    }
    done
fi
wait
echo "MATK done"