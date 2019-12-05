## Install m6APipe
There are two [Installation methods](https://github.com/kingzhuky/m6APipe/wiki/Installation) for m6APipe.

## Running test
### Conda
In the directory of m6APipe
```
nextflow main.nf -c nextflow.config --readPaths test_human --designfile designfile_test_single.csv -resume
```

## Running your own data
### Modify nextflow.config
#### General parameter
Edit the nextflow.config, modify genome file, analysis mode and aligners' index.
```
  fasta = "/home/zky/m6apipe/Genome/hg38/hg38_genome.fa"
  gtf = "/home/zky/m6apipe/Genome/hg38/hg38_genes.gtf"

  aligners = "star" // "star" OR "bwa" OR "tophat2" OR "hisat2" OR "none"
  peakCalling_mode = "independence" // "group" OR "independence"
  peakMerged_mode = "rank" // "rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc"
  expression_analysis_mode = "DESeq2" // "DESeq2" OR "EdgeR" OR "none"
  methylation_analysis_mode = "QNB" // "MATK" OR "QNB" OR "Wilcox-test" OR "MeTDiff"
```
#### Option parameter
Setting aligners' index according to aligner you chose( Default: "star" )
```
  tophat2_index = "/Path/to/Tophat2Index/*"
  hisat2_index = "/Path/to//Hisat2Index/*"
  bwa_index = "/Path/to/BWAIndex/*"
  star_index = "/Path/to/starindex"
```
### Setting designfile and comparefile
#### Designfile
Edit the nextflow.config and define "readPaths", "aligners", "designfile", "comparefile" and correspondiente alignment index for recommend.
Designfile is just like the following table with a comma (,) separated, which is .csv suffix file.

| Sample_ID | input_FileName | ip_FileName | Group |
| --- | --- | --- | --- |
| H1A_Endo | A | B | group_Endo |
| H1A_ES | C | D | group_ES |
| H1B_Endo | E | F | group_Endo |
| H1B_ES | G | H | group_ES |

>Tips
>1. A, B, C... mean the filenames of data, just like A.fastq.gz.
>2. If your data is .fastq.gz suffix file, please add the parameter of gzip, just like "--gzip true".
>3. If your filename of data is "Hela_cell_input.fastq.gz", please write its filename as "Hela_cell_input".

#### Comparefile
Comparefile is just like the following text which is a "\_vs\_" between two groups. 
>group_Endo_vs_group_ES

### Running command
#### data: fastq.gz
```
nextflow /path/to/m6APipe/main.nf -c /path/to/m6APipe/nextflow.config --readPaths /path/to/data/ --designfile /your/designfile.csv --comparefile /your/comparefile --gzip
```
#### data: fastq
```
nextflow /path/to/m6APipe/main.nf -c /path/to/m6APipe/nextflow.config --readPaths /path/to/data/ --designfile /your/designfile.csv --comparefile /your/comparefile --gzip false
```
#### data: bam
```
nextflow /path/to/m6APipe/main.nf -c /path/to/m6APipe/nextflow.config --readPaths /path/to/data/ --designfile /your/designfile.csv --comparefile /your/comparefile --aligners none
```