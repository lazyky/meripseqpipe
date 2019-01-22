#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : readsfilename
designfile=$1
reads_filename=$2
#Windows and linux newline ^M conversion
cat $designfile > tmp_designfile.txt 
dos2unix tmp_designfile.txt
#rename reads name by designfile
awk '{FS=",";OFS=""}{if (NR>1)print "mv ", $1, ".fastq ",   $1,"_",$2,"_",$3,"_",$4,".fastq"}' tmp_designfile.txt | grep ${reads_filename} | bash
#fastqc
mkdir fastqc
ls *.fastq | xargs -iLIST fastqc -o fastqc --noextract LIST