#!/bin/bash
## matk.sh tophat2 $matk_jar $designfile $gtf
## $1 argv 1 : uesd Aligner
## $2 argv 2 : matk_jar
## $3 argv 3 : designfile
## $4 argv 4 : gtf file
Aligner_name=$1
matk_jar=$2
designfile=$3
gtf_file=$4

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
do
    ls *input_${i}_*.bam | awk '{ORS=";"}{print $0}'
done
java -jar MATK-1.0.jar -peakCalling -ip "Mut1ip_treated_ip_1_tophat2_sort.bam;Mut2ip_treated_ip_1_tophat2_sort.bam" -input "Mut1input_treated_input_1_tophat2_sort.bam;Mut2input_treated_input_1_tophat2_sort.bam" -out test2.bed
java -jar MATK-1.0.jar -diff \
                           -control_ip "control_ip1.bam;control_ip2.bam;control_ip3.bam" \
                           -control_input "control_input1.bam;control_input2.bam;control_input3.bam" \
                           -treated_ip "treated_ip1.bam;treated_ip2.bam;treated_ip3.bam" \
                           -treated_input "treated_input1.bam;treated_input2.bam;treated_input3.bam" \
                           -control_bed control_peak.bed \
                           -treated_bed treated_peak.bed \
                           -gtf hg19.gtf \
                           -out m6A_differentiation.txt

java -jar MATK-1.0.jar -diff 
                           -control_ip "WT1ip_control_ip_0_tophat2_sort.bam;WT2ip_control_ip_0_tophat2_sort.bam" 
                           -control_input "WT1input_control_input_0_tophat2_sort.bam;WT2input_control_input_0_tophat2_sort.bam" 
                           -treated_ip "Mut1ip_treated_ip_1_tophat2_sort.bam;Mut2ip_treated_ip_1_tophat2_sort.bam" 
                           -treated_input "Mut1input_treated_input_1_tophat2_sort.bam;Mut2input_treated_input_1_tophat2_sort.bam" 
                           -control_bed test.bed 
                           -treated_bed test2.bed 
                           -gtf genes.gtf 
                           -out test.txt

PePr --chip1 Mut1ip_treated_ip_1_tophat2_sort.bam,Mut2ip_treated_ip_1_tophat2_sort.bam --input1 Mut1input_treated_input_1_tophat2_sort.bam,Mut2input_treated_input_1_tophat2_sort.bam \
     --chip2 WT1ip_control_ip_0_tophat2_sort.bam,WT2ip_control_ip_0_tophat2_sort.bam --input2 WT1input_control_input_0_tophat2_sort.bam,WT2input_control_input_0_tophat2_sort.bam \
     -f bam --diff