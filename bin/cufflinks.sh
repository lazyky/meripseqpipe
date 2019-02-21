#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : designfile
#$3 argv 3 : gtf file
#$4 argv 4 : THREAD_NUM
Aligner_name=$1
designfile=$2
gtf_file=$3
THREAD_NUM=$4
MAX_SITUATION=$(awk -F, '{if(NR>1)print int($4)}' $designfile | sort -r | head -1)
tag=$(for ((i=0;i<=$MAX_SITUATION;i++))
        do 
                if [ $i -lt $MAX_SITUATION ] 
                then echo -e "Situation"$i",\c"
                else echo -e "Situation"$i"\c"
                fi          
        done 
        )
bam_file_array=$(for ((i=0;i<=$MAX_SITUATION;i++));
        do 
               echo *input_${i}_${Aligner_name}*.bam | awk '{OFS=",";ORS=""}{for(x=1;x<NF;x++) print $x"," }END{print $x" "}'
        done 
        )
rm -f assembly_list_${Aligner_name}.txt
for bam_file in *input*${Aligner_name}*.bam
do
        cufflinks -p ${THREAD_NUM} -G ${gtf_file} --library-type fr-unstranded -o ${bam_file/.bam/} $bam_file
        echo "./"${bam_file/.bam/}"/transcripts.gtf" >> assembly_list_${Aligner_name}.txt
done
cuffmerge -o ./merged_gtf_${Aligner_name} -g ${gtf_file} -p ${THREAD_NUM} assembly_list_${Aligner_name}.txt
cuffdiff -o cuffdiff_${Aligner_name} \
         -L $tag \
         -p ${THREAD_NUM} \
         --time-series --multi-read-correct \
         --library-type fr-unstranded \
         ./merged_gtf_${Aligner_name}/merged.gtf ${bam_file_array}
