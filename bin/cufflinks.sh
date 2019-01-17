#!/bin/bash
#$1 argv 1 : uesd Aligner
awk -F, '{if( $2 == "control" && $3 == "input")print $1,$2,$3,$4}' OFS="_" ORS="_$1_sort.bam\n" tmp_designfile.txt | xargs samtools merge -f control_input_$1_sort.bam
awk -F, '{if( $2 == "control" && $3 == "ip")print $1,$2,$3,$4}' OFS="_" ORS="_$1_sort.bam\n" tmp_designfile.txt | xargs samtools merge -f control_ip_$1_sort.bam
awk -F, '{if( $2 == "treated" && $3 == "input")print $1,$2,$3,$4}' OFS="_" ORS="_$1_sort.bam\n" tmp_designfile.txt | xargs samtools merge -f treated_input_$1_sort.bam
awk -F, '{if( $2 == "treated" && $3 == "ip")print $1,$2,$3,$4}' OFS="_" ORS="_$1_sort.bam\n" tmp_designfile.txt | xargs samtools merge -f treated_ip_$1_sort.bam
macs2 callpeak -t control_ip_$1_sort.bam -c control_input_$1_sort.bam -g dm --outdir macs2_$1_control -n macs2_$1_control -p 1e-6 -f BAM --shift=150 --nomodel
macs2 callpeak -t treated_ip_$1_sort.bam -c treated_input_$1_sort.bam -g dm --outdir macs2_$1_treated -n macs2_$1_treated -p 1e-6 -f BAM --shift=150 --nomodel
mv macs2_$1_control/macs2_$1_control_summits.bed macs2_$1_control/macs2_$1_control.bed
mv macs2_$1_treated/macs2_$1_treated_summits.bed macs2_$1_treated/macs2_$1_treated.bed