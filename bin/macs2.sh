#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : designfile
Aligner_name=$1
designfile=$2
MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)

for ((i=0;i<=$MAX_SITUATION;i++))
do
    ls *input_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_situation_${i}_input_${Aligner_name}.bam
    ls *ip_${i}_${Aligner_name}*.bam | awk 'BEGIN{ORS=" "}{print $0}'|xargs samtools merge -f macs2_situation_${i}_ip_${Aligner_name}.bam
    macs2 callpeak -t macs2_situation_${i}_ip_${Aligner_name}.bam -c macs2_situation_${i}_input_${Aligner_name}.bam -g dm -n macs2_situation_${i}_${Aligner_name} -p 1e-6 -f BAM --shift=150 --nomodel
    mv macs2_situation_${i}_${Aligner_name}_summits.bed macs2_situation_${i}_${Aligner_name}.bed
done
