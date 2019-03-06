#!/bin/bash
## MATK_diffpeakCalling.sh tophat2 $matk_jar $designfile $gtf
## $1 argv 1 : TREATED_SITUATION_STARTPOINT
## $2 argv 2 : uesd Aligner
## $3 argv 3 : matk_jar
## $4 argv 4 : designfile
## $5 argv 5 : gtf file
TREATED_SITUATION_STARTPOINT=$1
Aligner_name=$2
matk_jar=$3
designfile=$4
gtf_file=$5

MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
for i in $(seq 0 $(expr ${TREATED_SITUATION_STARTPOINT} - 1));
do 
    for j in $(seq ${TREATED_SITUATION_STARTPOINT} ${MAX_SITUATION} )
    do
        control_ip_bam_file_array=$(echo *ip_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        control_input_bam_file_array=$(echo *input_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_ip_bam_file_array=$(echo *ip_${j}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        treated_input_bam_file_array=$(echo *input_${j}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
        java -jar ${matk_jar} -diff \
                        -control_ip "${control_ip_bam_file_array}" \
                        -control_input "${control_input_bam_file_array}" \
                        -treated_ip "${treated_ip_bam_file_array}" \
                        -treated_input "${treated_input_bam_file_array}" \
                        -control_bed MATK_peakCalling_situation_${i}.bed \
                        -treated_bed MATK_peakCalling_situation_${j}.bed \
                        -gtf ${gtf_file} \
                        -out MATK_diffpeakCalling_situation_${i}__${j}.txt
    done
done
