#!/bin/bash
## MATK_diffpeakCalling.sh tophat2 $matk_jar $designfile $gtf
## $1 argv 1 : matk_jar
## $2 argv 2 : designfile
## $3 argv 3 : gtf file
## $4 argv 4 : comparefile
matk_jar=$1
designfile=$2
gtf_file=$3
compare_str=$4
# Running MATK quantification
if [ compare_str != "two_group" ]; then
    # Running MATK quantification with compare_str
    group_id_1=$(echo $compare_str | awk 'BEGIN{FS="_vs_"}{print $1}')
    group_id_2=$(echo $compare_str | awk 'BEGIN{FS="_vs_"}{print $2}')
    echo $group_id_1    
    control_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    control_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_ip_bam_file_array=$(echo *ip_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_input_bam_file_array=$(echo *input_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    echo $compare_group
    java -jar ${matk_jar} -diff \
                    -control_ip "${control_ip_bam_file_array}" \
                    -control_input "${control_input_bam_file_array}" \
                    -treated_ip "${treated_ip_bam_file_array}" \
                    -treated_input "${treated_input_bam_file_array}" \
                    -control_bed *group_${group_id_1}.bed \
                    -treated_bed *group_${group_id_2}.bed \
                    -gtf ${gtf_file} \
                    -out MATK_diffm6A_${group_id_1}_${group_id_2}.txt
    echo $compare_group"  end "
else
    # Running MATK quantification without compare_str beacause of only two groups
    echo "no compare file"
    group_list=$(awk 'BEGIN{FS=","}NR>1{print $4}' $designfile |sort|uniq|awk 'BEGIN{ORS="\t"}{print $0}')
    group_id_1=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $1}')
    group_id_2=$(echo $group_list | awk 'BEGIN{FS="\t"}{print $2}')  
    echo $group_id_1   
    control_ip_bam_file_array=$(echo *ip_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    control_input_bam_file_array=$(echo *input_${group_id_1}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_ip_bam_file_array=$(echo *ip_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    treated_input_bam_file_array=$(echo *input_${group_id_2}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x";" }END{print $x""}')
    echo "group_list"
    java -jar ${matk_jar} -diff \
                    -control_ip "${control_ip_bam_file_array}" \
                    -control_input "${control_input_bam_file_array}" \
                    -treated_ip "${treated_ip_bam_file_array}" \
                    -treated_input "${treated_input_bam_file_array}" \
                    -control_bed *group_${group_id_1}.bed \
                    -treated_bed *group_${group_id_2}.bed \
                    -gtf ${gtf_file} \
                    -out MATK_diffm6A_${group_id_1}_${group_id_2}.txt
    echo $group_list"  end "
fi
wait
echo "diffMATK done"