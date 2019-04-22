#! /bin/sh
INPUTbam=$1
IPbam=$2
prefix=$3
chrName_file=$4
genomebin_dir=$5
THREAD_NUM=$6

#sort, index bam
#samtools sort $INPUTbam $INPUTbam.srt
#samtools sort $IPbam $IPbam.srt
#samtools index $INPUTbam.srt.bam
#samtools index $IPbam.srt.bam


#count total reads
samtools view -c $INPUTbam > $INPUTbam.readcounts.txt
samtools view -c $IPbam > $IPbam.readcounts.txt


echo "cal read counts for each bin"
mkdir $prefix.tmp
mkdir "$prefix.tmp/input"
awk -v bam="$INPUTbam" -v pre="$prefix" '{print "samtools view -b "bam" "$1 ">./"pre".tmp/input/"$1".bam"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD
awk -v pre="$prefix"  '{print "bamToBed -split -i < ./"pre".tmp/input/"$1".bam>./"pre".tmp/input/"$1".bed"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD
mkdir "$prefix.tmp/ip"
awk -v bam="$IPbam" -v pre="$prefix"  '{print "samtools view -b "bam" "$1 ">./"pre".tmp/ip/"$1".bam"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD
awk -v pre="$prefix"   '{print "bamToBed -split -i < ./"pre".tmp/ip/"$1".bam>./"pre".tmp/ip/"$1".bed"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD

awk -v pre="$prefix"  '{print "sortBed -i ./"pre".tmp/input/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/input/"$1".bin25.txt"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD
awk -v pre="$prefix"  '{print "sortBed -i ./"pre".tmp/ip/"$1".bed | intersectBed  -a '${genomebin_dir}'"$1".bin25.bed -b - -sorted -c > ./"pre".tmp/ip/"$1".bin25.txt"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD

echo "cal pval for each 25bp bin"
awk -v IPbam="$IPbam" -v INPUTbam="$INPUTbam" -v pre="$prefix" '{print "paste -d\\\"\\\\t\\\" ./"pre".tmp/ip/"$1".bin25.txt ./"pre".tmp/input/"$1".bin25.txt | perl cal_pval_for_bin.pl - "IPbam".readcounts.txt "INPUTbam".readcounts.txt> ./"pre".tmp/"$1".m6A.meyer.pval.txt"}' $chrName_file|xargs -iCMD -P$THREAD_NUM bash -c CMD

echo "bonferroni method to adjust the p value"
# ls $prefix.tmp/*.m6A.meyer.pval.txt|grep -v chrEBV|xargs -iF cat F > $prefix.human.m6A.meyer.pval.txt
# ls $prefix.tmp/*.m6A.meyer.pval.txt|grep chrEBV|xargs -iF cat F > $prefix.EBV.m6A.meyer.pval.txt
ls $prefix.tmp/*.m6A.meyer.pval.txt | xargs -iF cat F > $prefix.human.m6A.meyer.pval.txt

cat $prefix.human.m6A.meyer.pval.txt|wc -l|xargs -iCNT bash -c "awk '\$4<0.05/CNT{print \$0}' $prefix.human.m6A.meyer.pval.txt" > $prefix.human.m6A.meyer.IP.bonferroni_filter_0.05.txt
cat $prefix.human.m6A.meyer.pval.txt|wc -l|xargs -iCNT bash -c "awk '\$5<0.05/CNT{print \$0}' $prefix.human.m6A.meyer.pval.txt" > $prefix.human.m6A.meyer.INPUT.bonferroni_filter_0.05.txt

echo "cluster the bins and get the peaks: discard the concatenated windows that were smaller than 100 bp."

cmd1="perl cluster_bins.pl $prefix.human.m6A.meyer.IP.bonferroni_filter_0.05.txt |sortBed -i - > $prefix.human.m6A.meyer.IP.bonferroni_filter_0.05_peaks.bed"
cmd2="perl cluster_bins.pl $prefix.human.m6A.meyer.INPUT.bonferroni_filter_0.05.txt |sortBed -i - > $prefix.human.m6A.meyer.INPUT.bonferroni_filter_0.05_peaks.bed"

echo -e "$cmd1\n$cmd2"|xargs -iCMD -P$THREAD_NUM bash -c CMD
rm -r $prefix.tmp