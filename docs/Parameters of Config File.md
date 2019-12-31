The typical command for running the pipeline is as follows:
```
nextflow path/to/m6APipe/main.nf -c path/to/m6APipe/nextflow.config \
    --readPaths path/to/datapath \
    --designfile path/to/designfile \
    --comparefile path/to/comparefile \
    -profile standard[,docker] \
    -resume
```
As you can see, you can use the parameter like the above command "--fasta path/to/fasta_file" 
## Detail of Parameters

| Mandatory arguments|  |
| --- | --- |
| readPaths | Path to input data (must be surrounded with quotes) |
| fasta  | Path to genome sequence file ( .fa ) |
| gtf | Path to genome annotation file ( .gtf ) |
| designfile  | Path to file of designfile ( format: Sample_ID,input_FileName,ip_FileName,Group ) |
| comparefile  | Path to file of comparefile ( format: A_vs_B ) |
| aligners | star/bwa/tophat2/hisat2/none (must be surrounded with quotes) 'none' means input file is BAM file and skip alignment |
| peakCalling_mode  | "group" OR "independence" ( Group means that there are biological replicates in function of PeakCalling, while independence means no biological replicate) |
| peakMerged_mode |  "mspc" OR "rank" OR "macs2" OR "MATK" OR "metpeak" |
| methylation_analysis_mode | "MATK" OR "QNB" OR "Wilcox-test" OR "MeTDiff" OR "edgeR" OR "DESeq2" |
| expression_analysis_mode | "DESeq2" OR "edgeR" OR "none" |
| gzip | Boolean value ( true/false ) : True means your data is gzip compressed ( .gz ) while false is not |
| singleEnd  | Specifies that the input is single end reads |
| stranded | Specifies that the input is strand specific ( yes/no/reverse ), defalut is 'no' |
---
| Options parameters|  |
| --- | --- |
| tophat2_index | Path to tophat2 index, eg. "path/to/Tophat2Index/\*" |
| hisat2_index  | Path to hisat2 index, eg. "path/to/Hisat2Index/\*" |
| bwa_index  | Path to bwa index, eg. "path/to/BwaIndex/\*" |
| star_index | Path to star index, eg. "path/to/StarIndex/" |
| matk_jar | Path to the jar file of MATK, eg "path/to/MATK-1.0.jar" |
---
| Other parameters|  |
| --- | --- |
| outdir | The output directory where the results will be saved, defalut = $baseDir/results |
| email  | Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits |
| skip_qc | Skip all QC steps   |
| skip_expression | Skip all differential expression analysis steps |
| skip_peakCalling  | Skip all Peak Calling steps |
| skip_diffpeakCalling | Skip all Differential methylation analysis |
| skip_fastqc  | Skip FastQC |
| skip_fastp | Skip Fastp |
| skip_rseqc | Skip RSeQC   |
| skip_createbedgraph | Skip generating bigwig file for genebody coverage  |
| skip_metpeak  | Skip MeTPeak's process of Peak Calling steps |
| skip_macs2 | Skip MACS2's process of Peak Calling steps
| skip_matk  | Skip MATK's process of Peak Calling steps |
| skip_meyer | Skip Meyer's process of Peak Calling steps|
| skip_cufflinks | Skip the cufflinks process of differential expression analysis steps |
| skip_edger  | Skip the EdgeR process of differential expression analysis steps |
| skip_deseq2 | Skip the DESeq2 process of differential expression analysis steps |

### Designfile
Edit the nextflow.config and define "readPaths", "designfile", "comparefile", "aligners" and correspondiente alignment index for recommend.
Designfile is just like the following table with a comma (,) separated, which is .csv suffix file. You also can see in [designfile_test.csv]( https://github.com/kingzhuky/m6APipe/blob/master/test_data/designfile_test.csv).


| Sample_ID| input_FileName | ip_FileName |  Group |
| --- | --- | --- | --- |
| H1A_Endo | A | B | group_Endo |
| H1A_ES | C | D | group_ES |
| H1B_Endo | E | F | group_Endo |
| H1B_ES | G | H | group_ES |

>Tips
>1. A, B, C... mean the filenames of data, just like A.fastq.gz.
>2. If your data is .fastq.gz suffix file, please add the parameter of gzip, just like "--gzip true".
>3. If your filename of data is "Hela_cell_input.fastq.gz", please write its filename as "Hela_cell_input".
#### Sequence Data(.fastq)
While your data is paired-end, you need to make the file name format, just like "A_R1.fastq[.gz]" or  "A_1.fastq[.gz]".( A is the file name, like the above table's input_FileName ) 
For example, your paired-end data file name is the following text if your designfile is the above table.
```
# paired-end data file names
A_R1.fastq 
A_R2.fastq 
B_R1.fastq
B_R2.fastq
------  OR  -------
A_1.fastq 
A_2.fastq 
B_1.fastq
B_2.fastq
```

And then the single-end data file name is just like the following text.
```
# single-end data file names
A.fastq 
B.fastq
```
#### Aligned Data(.bam)
the bam data file name is just like the following text.
```
# single-end data file names
A.bam 
B.bam
```

### Comparefile
Comparefile is just like the following text which is a "\_vs\_" between two groups, just like the file [comparefile.txt](
https://github.com/kingzhuky/m6APipe/blob/master/test_data/comparefile.txt). 
>group_Endo_vs_group_ES

