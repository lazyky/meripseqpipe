#!/bin/bash
## MATK_quantification.sh $matk_jar $gtf $formatted_designfile ${task.cpus} ${peak_bed}
## $1 argv 1 : matk_jar
## $2 argv 2 : gtf file
## $3 argv 3 : designfile
## $4 argv 4 : merge_bed_file

matk_jar=$1
gtf_file=$2
designfile=$3
merge_bed_file=$4

sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
for sample_id in $sample_list
do
{
    ip_bam_file=$(ls ${sample_id}.ip*.bam)
    input_bam_file=$(ls ${sample_id}.input*.bam)
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