#!/bin/bash
## MATK_diffpeakCalling.sh tophat2 $matk_jar $designfile $gtf
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : comparefile
## $4 argv 4 : gtf file
## $5 argv 6 : THREAD_NUM
matk_jar=$1
designfile=$2
comparefile=$3
gtf_file=$4
THREAD_NUM=$5


#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done
if [ ${comparefile:=0} ]; then
    compare_list=$(cat $comparefile)
    for compare_group in $compare_list;
    do
    read -u 9
    {
        group_id_1=$(echo $compare_group | awk 'BEGIN{FS="_vs_"}{print $1}')
        group_id_2=$(echo $compare_group | awk 'BEGIN{FS="_vs_"}{print $2}')    
        control_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        control_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        java -jar ${matk_jar} -diff \
                        -control_ip "${control_ip_bam_file_array}" \
                        -control_input "${control_input_bam_file_array}" \
                        -treated_ip "${treated_ip_bam_file_array}" \
                        -treated_input "${treated_input_bam_file_array}" \
                        -control_bed MATK_group_${group_id_1}.bed \
                        -treated_bed MATK_group_${group_id_1}.bed \
                        -gtf ${gtf_file} \
                        -out diffMATK_group_${i}__${group_id_1}.bed
        echo >&9
    } &
    done
else
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS="\t"}{print $0}')
    group_id_1=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $1}')
    group_id_2=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $2}')  
    control_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    control_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    java -jar ${matk_jar} -diff \
                    -control_ip "${control_ip_bam_file_array}" \
                    -control_input "${control_input_bam_file_array}" \
                    -treated_ip "${treated_ip_bam_file_array}" \
                    -treated_input "${treated_input_bam_file_array}" \
                    -control_bed MATK_group_${group_id_1}.bed \
                    -treated_bed MATK_group_${group_id_1}.bed \
                    -gtf ${gtf_file} \
                    -out diffMATK_group_${i}__${group_id_1}.bed
fi
wait
echo "diffMATK done"