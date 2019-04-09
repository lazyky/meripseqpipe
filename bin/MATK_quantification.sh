#!/bin/bash
## MATK_quantification.sh $matk_jar $gtf $formatted_designfile ${task.cpus} ${peak_bed}
## $1 argv 1 : matk_jar
## $2 argv 2 : gtf file
## $3 argv 3 : designfile
## $4 argv 4 : THREAD_NUM
## $5 argv 5 : merge_bed_file

matk_jar=$1
gtf_file=$2
designfile=$3
THREAD_NUM=$4
merge_bed_file=$5

# if [ $flag_peakCallingbygroup -gt 0 ]; then
#     group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
#     for group_id in $group_list
#     do 
#     read -u 9
#     {
#         ip_bam_file_array=$(echo *ip_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
#         input_bam_file_array=$(echo *input_${group_id}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
#         java -jar $matk_jar -quantification \
#                            -ip "$ip_bam_file_array" \
#                            -input "$input_bam_file_array" \
#                            -bed $merge_bed_file \
#                            -gtf $gtf_file \
#                            -out MATK_group_${group_id}_quantification.bed
#         echo >&9
#     }&
#     done 
# else
sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for sample_id in $sample_list
do
{
    ip_bam_file=$(ls ${sample_id}.ip*.bam)
    input_bam_file=$(ls ${sample_id}.ip*.bam)
    java -jar $matk_jar -quantification \
                -ip "$ip_bam_file" \
                -input "$input_bam_file" \
                -bed $merge_bed_file \
                -gtf $gtf_file \
                -out MATK_${sample_id}_quantification.bed
    echo $sample_id > tmp.quantification.$sample_id
    awk 'BEGIN{FS="\t"}{print $5}' MATK_${sample_id}_quantification.bed >> tmp.quantification.$sample_id
}
done
echo "Peak_name" > tmp.MATK.quantification
awk 'BEGIN{FS="\t"}{print $4}' $merge_bed_file >> tmp.MATK.quantification
ls tmp.quantification.* |xargs paste tmp.MATK.quantification > MATK_quantification.matrix
wait
echo "MATK quantification done"