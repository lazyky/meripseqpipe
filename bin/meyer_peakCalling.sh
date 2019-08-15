#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : chrName.txt
#$3 argv 3 : genomebin_dir
#$4 argv 4 : peak_windows_number
#$4 argv 4 : THREAD_NUM
#$5 argv 5 : flag_peakCallingbygroup
designfile=$1
chrName_file=$2
genomebin_dir=$3
peak_windows_number=$4
THREAD_NUM=$5
flag_peakCallingbygroup=$6

function meyer_peakCalling()
{
    input_bam=$1
    ip_bam=$2
    prefix=$3
    chrName_file=$4
    genomebin_dir=$5
    peak_windows_number=$6
    THREAD_NUM=$7
    input_total_reads_count=$(samtools view -c $input_bam)
    ip_total_reads_count=$(samtools view -c $ip_bam)
    mkdir $prefix.tmp "$prefix.tmp/input" "$prefix.tmp/ip"
    awk -v bam="$input_bam" -v pre="$prefix" '
        {print "samtools view -b "bam" "$1 ">./"pre".tmp/input/"$1".bam; \
        bamToBed -split -i < ./"pre".tmp/input/"$1".bam>./"pre".tmp/input/"$1".bed; \
        sortBed -i ./"pre".tmp/input/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/input/"$1".bin25.txt"}' $chrName_file \
        | xargs -iCMD -P$THREAD_NUM bash -c CMD
    awk -v bam="$ip_bam" -v pre="$prefix"  '
        {print "samtools view -b "bam" "$1 ">./"pre".tmp/ip/"$1".bam; \
        bamToBed -split -i < ./"pre".tmp/ip/"$1".bam>./"pre".tmp/ip/"$1".bed; \
        sortBed -i ./"pre".tmp/ip/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/ip/"$1".bin25.txt"}' $chrName_file \
        | xargs -iCMD -P$THREAD_NUM bash -c CMD
    echo "cal pval for each 25bp bin"
    awk -v pre="$prefix" '
    {print "python meyer.py ./"pre".tmp/input/"$1".bin25.txt ./"pre".tmp/ip/"$1".bin25.txt '$input_total_reads_count' '$ip_total_reads_count' '$peak_windows_number' ./"pre".tmp/ip/"$1".m6A.meyer.pval.txt"}' $chrName_file \
    |xargs -iCMD -P$THREAD_NUM bash -c CMD
    cat $prefix.tmp/ip/*.m6A.meyer.pval.txt > meyer_${prefix}.bed
    awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$4}' meyer_${prefix}.bed > meyer_${prefix}_normalized.bed
    rm -rf $prefix.tmp
}

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
        samtools index ${group_id}_input.bam
        samtools index ${group_id}_ip.bam
        meyer_peakCalling ${group_id}_input.bam ${group_id}_ip.bam group_${group_id} $chrName_file $genomebin_dir $peak_windows_number $THREAD_NUM
    }
    done
else
    sample_list=$(awk 'BEGIN{FS=","}NR>1{print $1}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_id in $sample_list
    do
    {
        meyer_peakCalling ${sample_id}.input*.bam ${sample_id}.ip*.bam $sample_id $chrName_file $genomebin_dir $peak_windows_number $THREAD_NUM
    }
    done
fi
wait
echo "Meyer done"
