#!/bin/bash
## MATK_peakCalling.sh tophat2 $matk_jar $designfile
## $1 argv 1 : uesd Aligner
## $2 argv 2 : matk_jar
## $3 argv 3 : designfile
## $4 argv 4 : THREAD_NUM
Aligner_name=$1
matk_jar=$2
designfile=$3
THREAD_NUM=$4

#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
for ((i=1;i<=$MAX_SITUATION;i++))
do 
read -u 9
{
    ip_bam_file_array=$(echo *ip_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    input_bam_file_array=$(echo *input_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    java -jar $matk_jar -peakCalling -ip "$ip_bam_file_array" -input "$input_bam_file_array" -out MATK_situation_${i}_${Aligner_name}.bed
    echo >&9
}&
done 
