#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : designfile
#$2 argv 3 : THREAD_NUM
Aligner_name=$1
designfile=$2
THREAD_NUM=$3

#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
for ((i=0;i<=$MAX_SITUATION;i++))
do
read -u 9
{
    count=$(ls *input_${i}_${Aligner_name}*.bam| wc -w)
    if [ $count -gt 1 ]; 
    then
        ls *ip_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_situation_${i}_ip_${Aligner_name}.bam
        ls *input_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_situation_${i}_input_${Aligner_name}.bam
    else 
        ls *ip_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print "mv "$0," macs2_situation_'${i}'_ip_'${Aligner_name}'.bam"}' | bash
        ls *input_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print "mv "$0," macs2_situation_'${i}'_input_'${Aligner_name}'.bam"}' | bash
    fi
    macs2 callpeak -t macs2_situation_${i}_ip_${Aligner_name}.bam -c macs2_situation_${i}_input_${Aligner_name}.bam -g hs -n macs2_situation_${i}_${Aligner_name} -p 1e-6 -f BAM --nomodel
    mv macs2_situation_${i}_${Aligner_name}_summits.bed macs2_situation_${i}_${Aligner_name}.bed
    echo >&9
}&
done
