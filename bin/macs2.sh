#!/bin/bash
#bash macs2.sh <designfile> <THREAD_NUM> <whether_peakcallingbygroup>
#$1 argv 1 : designfile
#$2 argv 2 : genome_size
#$3 argv 3 : flag_peakCallingbygroup
#$4 argv 4 : THREAD_NUM
## flag_peakCallingbygroup: 1(group) 0(sample)
designfile=$1
genome_size=$2
flag_peakCallingbygroup=$3
THREAD_NUM=$4
arguments=`echo ${@:5}`
# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9  
done

# Check the mode of peakcalling
if [ $flag_peakCallingbygroup -gt 0 ]; then
    ## Running Macs2 for peakcalling per group
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        #check whether the number of sample in the group > 2
        count=$(ls *input_${group_id}*.bam| wc -w)
        if [ $count -gt 1 ]; then
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_ip.bam
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_input.bam
        else 
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_ip.bam"}' | bash
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_input.bam"}' | bash
        fi
        macs2 callpeak -t macs2_group_${group_id}_ip.bam -c macs2_group_${group_id}_input.bam -g $genome_size -n macs2_group_${group_id} $arguments -f BAM --nomodel
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"10^-$8}' macs2_group_${group_id}_peaks.narrowPeak > macs2_group_${group_id}_normalized.bed
        mv macs2_group_${group_id}_summits.bed macs2_group_${group_id}.summits
        echo >&9
    }&
    done
else
    ## Running Macs2 for peakcalling per samole
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    read -u 9
    {
        macs2 callpeak -t ${sample_id}.ip*.bam -c ${sample_id}.input*.bam -g $genome_size -n macs2_${sample_id} $arguments -f BAM --nomodel
        awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"10^-$8}' macs2_${sample_id}_peaks.narrowPeak > macs2_${sample_id}_normalized.bed
        mv macs2_${sample_id}_summits.bed macs2_${sample_id}.summits
        echo >&9
    }&
    done
fi
wait
echo "Macs2 done"