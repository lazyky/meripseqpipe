#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3

#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done

if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        count=$(ls *input_${group_id}*.bam| wc -w)
        if [ $count -gt 1 ]; then
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_ip.bam
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_group_${group_id}_input.bam
        else 
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_ip.bam"}' | bash
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," macs2_group_'${group_id}'_input.bam"}' | bash
        fi
        macs2 callpeak -t macs2_group_${group_id}_ip.bam -c macs2_group_${group_id}_input.bam -g hs -n macs2_group_${group_id} -p 1e-6 -f BAM --nomodel
        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' macs2_group_${group_id}_peaks.narrowPeak > macs2_group_${group_id}.bed
        mv macs2_group_${group_id}_summits.bed macs2_group_${group_id}.summits
        echo >&9
    }&
    done
else
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    read -u 9
    {
        macs2 callpeak -t ${sample_id}.ip*.bam -c ${sample_id}.input*.bam -g hs -n macs2_${sample_id} -p 1e-6 -f BAM --nomodel
        awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' macs2_${sample_id}_peaks.narrowPeak > macs2_${sample_id}.bed
        mv macs2_${sample_id}_summits.bed macs2_${sample_id}.summits
        echo >&9
    }&
    done
fi
wait
echo "Macs2 done"