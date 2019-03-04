#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : filename
designfile=$1
filename=$2
filesuffix=$3
#Windows and linux newline ^M conversion
cat $designfile > tmp_designfile.txt 
dos2unix tmp_designfile.txt
#rename reads name by designfile
if [ $filesuffix == "bam" ]; then 
    awk '{FS=",";OFS=" "}{if (NR>1)print "mv "$1".bam",$1"_"$2"_"$3"_"$4"_aligners.bam"}' tmp_designfile.txt | uniq | grep $filename | bash
elif [ $filesuffix == "fastq" ]; then
    awk '{FS=",";OFS=" "}{if (NR>1)print "mv "$1".fastq",$1"_"$2"_"$3"_"$4"_aligners.fastq"}' tmp_designfile.txt | uniq | grep $filename | bash
else
    echo "you had enter unexpected word"
fi
