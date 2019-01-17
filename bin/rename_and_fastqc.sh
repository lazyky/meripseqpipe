#!/bin/bash
#$1 argv 1 : designfile
#$2 argv 2 : readsfilename
cat $1 > tmp_designfile.txt 
dos2unix tmp_designfile.txt 
awk '{FS=",";OFS=""}{if (NR>1)print "mv ",$1,".fastq ",$1,"_",$2,"_",$3,"_",$4,".fastq"}' tmp_designfile.txt | grep $2 | bash
mkdir fastqc
ls *.fastq | xargs -iLIST fastqc -o fastqc --noextract LIST