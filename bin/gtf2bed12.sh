#!/bin/bash
# gtfToGenePred were download form UCSC  : http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
#$1 argv 1 : gtf file, the name of output is genes.bed
gtfToGenePred -genePredExt -geneNameAsName2 $1 genes.tmp
awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' genes.tmp >  ${1/.gtf/.bed}
rm genes.tmp
