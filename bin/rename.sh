#!/bin/bash
#$1 argv 1 : Aligners name
#$2 argv 2 : designfile
Aligners_name=$1
designfile=$2
#rename reads name by designfile
if [ $Aligners_name == "aligners" ]; then 
    awk '{FS=",";OFS=" "}{if (NR>1)print "ln "$1".bam",$1"_"$2"_"$3"_"$4"_'${Aligners_name}'.bam"}' $designfile| bash
else
    awk '{FS=",";OFS=" "}{if (NR>1)print "ln "$1"_aligners_'${Aligners_name}'.bam",$1"_"$2"_"$3"_"$4"_'${Aligners_name}'.bam"}' $designfile | bash
fi
wait
    echo "File name is renamed by designfile"