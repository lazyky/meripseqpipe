#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3

# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9 
done

function sort_and_transferbed()
{
    bed_file=$1
    bed_anno_file=$2
    outdir=$3
    sort -k5,5 -n -r ${bed_file} | awk '{ print $1":"$2"-"$3}' > ${outdir}/tmp.${bed_file}
    cat ${outdir}/tmp.${bed_file} | xargs -iID grep ID ${bed_anno_file} | awk '{print $2}' > ${outdir}/tmp.${bed_file}.location
    #rm -rf tmp.${outdir}
}
function mergebed_by_rank()
{
    prefix_id=$1
    out_prefix=$2
    mkdir tmp.${out_prefix}
    cat *${prefix_id}*.bed | awk '{print $1"\t"$2*1"\t"$3*1"\t"$1":"$2"-"$3}' > tmp.${out_prefix}/bedtools_${prefix_id}_all_peaks
    sortBed -i tmp.${out_prefix}/bedtools_${prefix_id}_all_peaks |mergeBed -i - -c 4,4 -o collapse,count | awk '{print $1"\t"$2"\t"$3"\t"$1":"$2"-"$3"\t"$4}'  > tmp.${out_prefix}/bedtools_${prefix_id}
    awk -F "\t" '{print $4,$5}' tmp.${out_prefix}/bedtools_${prefix_id} | awk -F '[," "]+' '{for (i=2 ;i<=NF;i++) printf $i" "$1"\n" }' | sort -k1 | uniq > tmp.${out_prefix}/bed_anno_file
    for bedfile in *${prefix_id}*.bed
    do
        sort_and_transferbed $bedfile tmp.${out_prefix}/bed_anno_file tmp.${out_prefix}
    done
    paste -d "\t" tmp.${out_prefix}/tmp*location > ${out_prefix}.bedlist
    peak_number=$(wc -l tmp.${out_prefix}/bed_anno_file | cut -d " " -f 1)
    Rscript merge_peaks_by_rank.R ${out_prefix}.bedlist ${peak_number} ${out_prefix}.bed
}

# if the number of peakcalling tools > 2
if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        mergebed_by_rank ${group_id} merged_group_${group_id}
        echo >&9
    }&
    done
else
    sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for sample_group_id in ${sampleinfo_list}
    do
    read -u 9
    {
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        ## Adding the information of group
        for samplefile in *${sample_id}*.bed
        do
            mv $samplefile ${samplefile/_normalized.bed/}_${group_id}_normalized.bed
        done
        mergebed_by_rank ${sample_id} merged_${sample_id}
        echo >&9
    }&
    done
    wait
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        mergebed_by_rank ${group_id} merged_group_${group_id}
        echo >&9
    }&
    done
fi
wait
mergebed_by_rank normalized Rankmerged_peaks
echo "bedtools merged peaks done"


