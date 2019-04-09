#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : THREAD_NUM
#$3 argv 3 : merge_bed_file
#$4 argv 4 : output_bam_stat_file
designfile=$1
THREAD_NUM=$2
merge_bed_file=$3
output_bam_stat_file=$4

mkfifo tmp
exec 9<>tmp
#rm -rf /tmp

for ((i=1;i<=${THREAD_NUM:=1};i++))
do
    echo >&9
done
sampleinfo_list=$(awk 'BEGIN{FS=","}NR>1{print $1","$4}' $designfile |sort|uniq|awk 'BEGIN{ORS=" "}{print $0}')
echo "Total_Reads" > $output_bam_stat_file
for sample_group_id in ${sampleinfo_list}
do
read -u 9
{
    sample_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $1}')
    group_id=$(echo ${sample_group_id} | awk 'BEGIN{FS=","}{print $2}')
    
    input_bam_file=$(ls ${sample_id}.input*.bam | awk '{ORS=" "}{print $0}')
    echo -e ${input_bam_file}"\t" | awk 'BEGIN{ORS=""}{print $0}' > ${sample_id}.bam_stat.txt
    samtools view -c ${input_bam_file} >> ${sample_id}.bam_stat.txt

    #Setting colnames of peaks input count
    echo ${sample_id}.input*.bam \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\tPeakName\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${merge_bed_file}.${group_id}.${sample_id}.input.count
    #Count input peaks
    bedtools multicov -bams ${input_bam_file} -bed ${merge_bed_file} >> ${merge_bed_file}.${group_id}.${sample_id}.input.count
    
    ip_bam_file=$(ls ${sample_id}.ip*.bam | awk '{ORS=" "}{print $0}')
    echo -e ${ip_bam_file}"\t" | awk 'BEGIN{ORS=""}{print $0}' >> ${sample_id}.bam_stat.txt
    samtools view -c ${ip_bam_file} >> ${sample_id}.bam_stat.txt
    #Setting colnames of peaks ip count
    echo ${sample_id}.ip*.bam \
    | awk 'BEGIN{ORS=""}{print "chrom\tchromStart\tchromEND\tPeakName\t"}{for(x=1;x<NF;x++) print $x"\t" }END{print $x"\n"}' \
    > ${merge_bed_file}.${group_id}.${sample_id}.ip.count
    #Count ip peaks
    bedtools multicov -bams ${ip_bam_file} -bed ${merge_bed_file} >> ${merge_bed_file}.${group_id}.${sample_id}.ip.count
}&
done
wait
cat *.bam_stat.txt >> $output_bam_stat_file
echo "bedtools count done"