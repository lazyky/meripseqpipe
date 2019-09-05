#!/bin/bash
#$1 argv 1 : peakCalling_tools_count
peakCalling_tools_main=$1
peakCalling_tools_count=$2

cat ${peakCalling_tools_main}*normalized.bed | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean > ${peakCalling_tools_main}_allPeaks.bed
ls *normalized.bed |grep -v ${peakCalling_tools_main} | xargs -i cat {} | sortBed -i - | mergeBed -i - -c 4,5 -o count,mean > others_allPeaks.bed
intersectBed -a ${peakCalling_tools_main}_allPeaks.bed -b others_allPeaks.bed -u | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$5}' > bedtools_merged_allpeaks.bed