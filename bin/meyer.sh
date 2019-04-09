#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3

if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    {
        count=$(ls *input_${group_id}*.bam| wc -w)
        if [ $count -gt 1 ]; then
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_ip.bam
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_input.bam
        else 
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_ip.bam"}' | bash
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_input.bam"}' | bash
        fi
        bash meyer_m6A_peak_call_by_yeying.sh macs2_group_${group_id}_input.bam macs2_group_${group_id}_ip.bam macs2_group_${group_id} $THREAD_NUM
    }
    done
else
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    {
        bash meyer_m6A_peak_call_by_yeying.sh ${sample_id}.input*.bam  ${sample_id}.ip*.bam macs2_${sample_id} $THREAD_NUM
    }
    done
fi
wait
echo "Meyer done"