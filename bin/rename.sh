#!/bin/bash
#bash rename.sh 

#$1 argv 1 : Aligners name
#$2 argv 2 : designfile
Aligners_name=$1
designfile=$2
sampleinfo_list=$(awk 'BEGIN{FS=",";OFS=","}NR>1{print $1,$2,$3,$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
# Rename the name of bamfiles
for sample_group_id in ${sampleinfo_list}
do
{
    sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
    group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $4}')
    input_sample_name=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
    ip_sample_name=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $3}')
    if [ $Aligners_name == "none" ]; then 
        ln -f ${input_sample_name}".bam" ${sample_id}".input_"${group_id}".bam"
        if [ $? != "0" ]; then
            echo "You may check your designfile, because there is something wrong in the process of rename process"
            exit 1
        fi
        ln -f ${ip_sample_name}".bam" ${sample_id}".ip_"${group_id}".bam"
        if [ $? != "0" ]; then
            echo "You may check your designfile, because there is something wrong in the process of rename process"
            exit 1
        fi
    else
        ln -f ${input_sample_name}"_"$Aligners_name".bam" ${sample_id}".input_"${group_id}".bam"
        if [ $? != "0" ]; then
            echo "You may check your designfile, because there is something wrong in the process of rename process"
            exit 1
        fi
        ln -f ${ip_sample_name}"_"$Aligners_name".bam" ${sample_id}".ip_"${group_id}".bam"
        if [ $? != "0" ]; then
            echo "You may check your designfile, because there is something wrong in the process of rename process"
            exit 1
        fi
    fi
}
done
wait
echo "File name is renamed by designfile"