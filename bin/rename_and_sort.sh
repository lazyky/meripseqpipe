#!/bin/bash
#$1 argv 1 : uesd Aligner
awk '{FS=",";OFS=""}{if (NR>1)print "mv ",$1,"_'$1'.bam ",$1,"_",$2,"_",$3,"_",$4,"_'$1'.bam"}' tmp_designfile.txt | bash
awk '{FS=",";OFS=""}{if (NR>1)print "samtools sort ",$1,"_",$2,"_",$3,"_",$4,"_'$1'.bam ","-o ",$1,"_",$2,"_",$3,"_",$4,"_'$1'_sort.bam"}' tmp_designfile.txt | bash
ls *$1_sort.bam | xargs -I@ samtools index @ 