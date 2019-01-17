#!/bin/bash
#$1 argv 1 : uesd Aligner
#$2 argv 2 : bed12
for bam_rseqc in *$1*.bam
do 
    infer_experiment.py -i $bam_rseqc -r $2 > ${bam_rseqc/_sort.bam/}.infer_experiment.txt
    junction_annotation.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r $2
    bam_stat.py -i $bam_rseqc > ${bam_rseqc/_sort.bam/}.bam_stat.txt
    junction_saturation.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r $2 2> ${bam_rseqc/_sort.bam/}.junction_annotation_log.txt
    inner_distance.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.rseqc -r $2
    read_distribution.py -i $bam_rseqc -r $2 > ${bam_rseqc/_sort.bam/}.read_distribution.txt
    read_duplication.py -i $bam_rseqc -o ${bam_rseqc/_sort.bam/}.read_duplication
done
