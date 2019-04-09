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
# if the number of peakcalling tools > 2
if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        cat *${group_id}*.bed | awk '{print $1"\t"$2*1"\t"$3*1"\t"$1":"$2"-"$3}' > bedtools_group_${group_id}_all_peaks #"_group_" for devive
        sortBed -i bedtools_group_${group_id}_all_peaks |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > bedtools_group_${group_id}.bed
        echo >&9
    }&
    done
else
    sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    # if the number of peakcalling tools > 2
    for sample_group_id in ${sampleinfo_list}
    do
    read -u 9
    {
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        cat *${sample_id}*.bed | awk '{print $1"\t"$2*1"\t"$3*1"\t"$1":"$2"-"$3}' > bedtools_${sample_id}_all_peaks
        sortBed -i bedtools_${sample_id}_all_peaks |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > bedtools_${group_id}_${sample_id}.bed
        echo >&9
    }&
    done
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        cat *${group_id}*.bed | awk '{print $1"\t"$2*1"\t"$3*1"\t"$1":"$2"-"$3}' > bedtools_group_${group_id}_all_peaks #"_group_" for devive
        sortBed -i bedtools_group_${group_id}_all_peaks |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > bedtools_group_${group_id}.bed
        echo >&9
    }&
    done
fi
cat bedtools_group*all_peaks | sortBed -i - |mergeBed -i - -c 4,4 -o collapse,count | awk '$5>1{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3}' > bedtools_merged_peaks.bed
wait
echo "bedtools merged peaks done"