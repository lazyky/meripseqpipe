#!/bin/bash
#$1 argv 1 : aligners_tools_name
#$2 argv 2 : peak_calling_tools_name
#$3 argv 3 : designfile
#$4 argv 4 : THREAD_NUM
aligners_tools_name=$1
peak_calling_tools_name=$2
designfile=$3
THREAD_NUM=$4
mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
peak_calling_merged_bed_file=${peak_calling_tools_name}_merged_peaks.bed
for ((i=1;i<=$MAX_SITUATION;i++))
do 
read -u 9
{   
    input_bam_file_array=$(ls *input_${i}_${aligners_tools_name}*.bam | awk '{ORS=" "}{print $0}')
    #Setting colnames of peaks input count
    echo *input_${i}_${aligners_tools_name}*.bam \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\tPeakName\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${aligners_tools_name}_${peak_calling_tools_name}_situation_${i}_input.count
    #Count input peaks
    bedtools multicov -bams ${input_bam_file_array} -bed ${peak_calling_merged_bed_file} >> ${aligners_tools_name}_${peak_calling_tools_name}_situation_${i}_input.count
    
    ip_bam_file_array=$(ls *ip_${i}_${aligners_tools_name}*.bam | awk '{ORS=" "}{print $0}')
    #Setting colnames of peaks ip count
    echo *ip_${i}_${aligners_tools_name}*.bam \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${aligners_tools_name}_${peak_calling_tools_name}_situation_${i}_ip.count
    #Count ip peaks
    bedtools multicov -bams ${ip_bam_file_array} -bed ${peak_calling_merged_bed_file} >> ${aligners_tools_name}_${peak_calling_tools_name}_situation_${i}_ip.count
}&
done 

wait
echo "done"
exec 9<&-
exec 9>&-
