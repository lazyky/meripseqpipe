#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
#$4 argv 4 : peakCalling_tools_count
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3
peakCalling_tools_count=$4

# Define a multi-threaded run channel
mkfifo tmp
exec 9<>tmp
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9 
done

function mergebed_for_BioRepeat()
{
    bed_file=$1
    bed_anno_file=$2
    outdir=$3
    dotnet CLI.dll -i rep1.bed -i rep2.bed -r bio -s 1E-8 -w 1E-4 -o 
}
function mergebed_for_TecRepeat()
{
    prefix_id=$1
    out_prefix=$2
    peakCalling_tools_count=$3
    mkdir tmp.${out_prefix}
    bedfile_array=$(ls *prefix_id*.bed | awk '{ORS=" "}{print "-i",$0}')
    dotnet CLI.dll $bedfile_array -r tec -c $peakCalling_tools_count -s 1E-8 -w 1E-4 -o Tec_$prefix_id
    mv Tec_$prefix_id

}
    ln -s !{mspc_dir} ./
    ls exomePeak*.bed | awk '{ORS=" "}{print "-i",$0}'| awk '{print "dotnet CLI.dll",$0,"-r bio -w 1E-4 -s 1E-8"}' | bash
    ls metpeak*.bed | awk '{ORS=" "}{print "-i",$0}'| awk '{print "dotnet CLI.dll",$0,"-r bio -w 1E-4 -s 1E-8"}' | bash
    for bed in */*.bed
    do
        mv $bed ${bed/%ConsensusPeaks.bed/temp.bed}
    done
    ls */*.bed | awk '{ORS=" "}{print "-i",$0}'| awk '{print "dotnet CLI.dll",$0,"-r bio -w 1E-4 -s 1E-8"}' | bash
# if the number of peakcalling tools > 2
if [ $flag_peakCallingbygroup -gt 0 ]; then
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebed_by_rank ${group_id} merged_group_${group_id}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${group_id}*.bed |sortBed -i - > merged_group_${group_id}.bed
        fi
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
        if [ $peakCalling_tools_count -gt 1 ]; then
            mergebed_by_rank ${sample_id} merged_${sample_id}
        else
            awk '{OFS="\t";$5=10^-$5;print }' *${sample_id}*.bed |sortBed -i - > merged_${sample_id}.bed
        fi
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
mergebed_by_rank merged_group Rankmerged_peaks
echo "bedtools merged peaks done"


