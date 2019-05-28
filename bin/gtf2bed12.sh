#!/bin/bash
#$1 argv 1 : gtf file, the name of output is genes.bed
gtf_file=$1
gtfToGenePred -genePredExt -geneNameAsName2 ${gtf_file} genes.tmp
genePredToBed genes.tmp ${gtf_file/.gtf/.bed}
rm genes.tmp