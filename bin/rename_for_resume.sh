#!/bin/bash
#bash rename.sh 
#$1 argv 1 : designfile
designfile=$1
sampleinfo_list=$(awk 'BEGIN{FS=",";OFS=","}NR>1{print $1,$2,$3,$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
# Rename the name of bamfiles
mkdir bam_for_resume
for sample_group_id in ${sampleinfo_list}
do
{
    sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
    group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $4}')
    input_sample_name=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
    ip_sample_name=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $3}')
    ln -f ${sample_id}".input_"${group_id}"_sort.bam" bam_for_resume/${input_sample_name}".bam"
    ln -f ${sample_id}".ip_"${group_id}"_sort.bam" bam_for_resume/${ip_sample_name}".bam"
}
done
wait
echo "File name is renamed by designfile"