#!/bin/bash
## MATK_diffpeakCalling.sh tophat2 $matk_jar $designfile $gtf
## $1 argv 1 : TREATED_SITUATION_STARTPOINT
## $2 argv 2 : uesd Aligner
## $3 argv 3 : matk_jar
## $4 argv 4 : designfile
## $5 argv 5 : gtf file
## $6 argv 6 : THREAD_NUM
TREATED_SITUATION_STARTPOINT=$1
Aligner_name=$2
matk_jar=$3
designfile=$4
gtf_file=$5
THREAD_NUM=$6

#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
for i in $(seq 1 $(expr ${TREATED_SITUATION_STARTPOINT} - 1));
do
    for j in $(seq ${TREATED_SITUATION_STARTPOINT} ${MAX_SITUATION} )
    do
    read -u 9
    {
        control_ip_bam_file_array=$(echo *ip_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        control_input_bam_file_array=$(echo *input_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_ip_bam_file_array=$(echo *ip_${j}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_input_bam_file_array=$(echo *input_${j}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        java -jar ${matk_jar} -diff \
                        -control_ip "${control_ip_bam_file_array}" \
                        -control_input "${control_input_bam_file_array}" \
                        -treated_ip "${treated_ip_bam_file_array}" \
                        -treated_input "${treated_input_bam_file_array}" \
                        -control_bed MATK_situation_${i}_${Aligner_name}.bed \
                        -treated_bed MATK_situation_${j}_${Aligner_name}.bed \
                        -gtf ${gtf_file} \
                        -out diffMATK_situation_${i}__${j}_${Aligner_name}.bed
        echo >&9
    } &
   done
done
