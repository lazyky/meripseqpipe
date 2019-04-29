#$/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : flag_peakCallingbygroup
designfile=$1
THREAD_NUM=$2
flag_peakCallingbygroup=$3

#定义描述符为9的管道
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp                   #关联后的文件描述符拥有管道文件的所有特性,所以这时候管道文件可以删除，我们留下文件描述符来用就可以了
for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9                   #&9代表引用文件描述符9，这条命令代表往管道里面放入了一个"令牌"
done
function mspc_merged_peaks_by_tagID()
{
    tag_id=$1
    prefix=$2
    mode=$3
    count=$(ls *${tag_id}*.bed | wc -w)
    if [ $count -gt 1 ]; then
        ls *${tag_id}*.bed |awk '{ORS=" "}{print "-i",$0}' | awk '{print "dotnet CLI.dll",$0,"-r '${mode}' -w 1E-4 -s 1E-8 > '${tag_id}'.log "}' | bash
        mv */ConsensusPeaks.bed ${prefix}
    else 
        ln *${tag_id}*.bed ${prefix}
    fi
}

# if the number of peakcalling tools > 2
if [ $flag_peakCallingbygroup -gt 0 ]; then
    ln -s ${mspc_directory}/* ./
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    for group_id in $group_list
    do
    read -u 9
    {
        mspc_merged_peaks_by_tagID group_${group_id} mspc_group_${group_id}.bed Tec
        echo >&9
    }&
    done
else
    sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
    # if the number of peakcalling tools > 2
    for sample_group_id in ${sampleinfo_list}
    do
    read -u 9
    {
        sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        mspc_merged_peaks_by_tagID ${sample_id} ${group_id}_${sample_id}.bed Tec 
        echo >&9
    }&
    done
    for sample_group_id in ${sampleinfo_list}
    do
    read -u 9
    {
        group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
        mspc_merged_peaks_by_tagID group_${group_id} mspc_group_${group_id}.bed Bio
        echo >&9
    }&
    done
fi
wait
mspc_merged_peaks_by_tagID mspc_group_${group_id} mspc_merged_peaks.bed Bio
echo "mspc merged peaks done"
 