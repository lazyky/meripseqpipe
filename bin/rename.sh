#!/bin/bash
#$1 argv 1 : Aligners name
#$2 argv 2 : designfile
Aligners_name=$1
designfile=$2
#rename reads name by designfile
if [ $Aligners_name == "none" ]; then 
    awk '{FS=",";OFS=" "}{if (NR>1)print "ln "$2".bam",$1".input_"$4".bam; ln "$3".bam",$1".ip_"$4".bam"}' $designfile | bash
else
    awk '{FS=",";OFS=" "}{if (NR>1)print "ln "$2"_'$Aligners_name'.bam",$1".input_"$4".bam; ln "$3"_'$Aligners_name'.bam",$1".ip_"$4".bam"}' $designfile | bash
fi
wait
echo "File name is renamed by designfile"