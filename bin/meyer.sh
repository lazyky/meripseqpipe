#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : chrName.txt
#$3 argv 3 : genomebin_dir
#$4 argv 4 : THREAD_NUM
#$5 argv 5 : flag_peakCallingbygroup
designfile=$1
chrName_file=$2
genomebin_dir=$3
THREAD_NUM=$4
flag_peakCallingbygroup=$5
#head -24 /data1/Database/hg38/hg38.chrom.sizes >hg38.chrom25.sizes
#grep -w 'chrM' /data1/Database/hg38/hg38.chrom.sizes >>hg38.chrom25.sizes
#sed -i 's/chrM/chrMT/' hg38.chrom25.sizes
#genomeLengthFile="hg38.chrom25.sizes"

#bedtools makewindows -g $genomeLengthFile -w 25 > genome.bin25.bed
#bedtools sort -i genome.bin25.bed >genome.bin25.srt.bed
#wc -l genome.bin25.bed | cut -d " " -f 1 >genome.bin25.count


if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    {
        count=$(ls *input_${group_id}*.bam| wc -w)
        if [ $count -gt 1 ]; then
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f ${group_id}_ip.bam
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f ${group_id}_input.bam
        else 
            ls *ip_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," '${group_id}'_ip.bam"}' | bash
            ls *input_${group_id}*.bam | awk 'BEGIN{ORS=" "}{print "ln "$0," '${group_id}'_input.bam"}' | bash
        fi
        bash meyer_m6A_peak_call_by_yeying.sh ${group_id}_input.bam ${group_id}_ip.bam $group_id $chrName_file $genomebin_dir $THREAD_NUM
    }
    done
else
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    {
        bash meyer_m6A_peak_call_by_yeying.sh ${sample_id}.input*.bam ${sample_id}.ip*.bam $sample_id $THREAD_NUM
    }
    done
fi
wait
echo "Meyer done"
