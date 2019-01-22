#!/bin/bash
#$1 argv 1 : gtf file, the name of output is genes.bed
gtf_file=$1
gtfToGenePred -genePredExt -geneNameAsName2 ${gtf_file} genes.tmp
awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' genes.tmp >  ${gtf_file/.gtf/.bed}
rm genes.tmp
