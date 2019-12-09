
## m6APipe
**MeRIP-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)
### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker/singularity containers making installation trivial and results highly reproducible. N6-methyladenosine (m6A) is the most prevalent modification in the mRNA of many eukaryotic species, including yeast, plants, flies, and mammals. In order to analyze m6A-seq data, we developed a user-friendly, integrated analysis pipeline called m6APipe based on Nextflow. It integrated ten main functional modules including preprocessing, QC, read mapping, peak calling, merging peaks, differential methylation analysis, differential expression analysis, motif search, annotation, and data visualization. 
### Documentation   
A full tutorial of m6APipe can be found at Wiki page of this project. plz go to the https://github.com/kingzhuky/m6APipe/wiki
#### Quickstart
##### Install Nextflow
Install Nextflow by using the following command:
```
curl -s https://get.nextflow.io | bash 
```
##### Install m6APipe
```
git clone https://github.com/kingzhuky/m6APipe.git
```
##### Build environment
Building environment by docker
```
docker pull kingzhuky/m6apipe
```
Or Building environment by conda. See details in [Installation](https://github.com/kingzhuky/m6APipe/wiki/Installation).
##### Launch m6APipe
```
nextflow run path/to/m6APipe/main.nf -c nextflow.config -profile docker
```
### Pipeline Description

#### Input files

The m6APipe pipeline needs as the input following files:
* Paired INPUT and IP samples reads, *.fastq or *.fastq.gz
* Genome assembly, *.fa
* Genome annotation, *.gtf
* Designfile contains grouping and pairing information , *.csv ( comma separated, "," ) 
* Comparefile for differential analysis, *.txt ( "_vs_" separeted, "," ) 


##### Designfile
Edit the nextflow.config and define "readPaths", "designfile", "comparefile", "aligners" and correspondiente alignment index for recommend.
Designfile is just like the following table with a comma (,) separated, which is .csv suffix file. You also can see in [designfile_test.csv]( https://github.com/kingzhuky/m6APipe/blob/master/test_data/designfile_test.csv).

| Sample_ID| input_FileName |ip_FileName |  Group |
| - | -| - | - |
| H1A_Endo | A | B | group_Endo |
| H1A_ES | C | D | group_ES |
| H1B_Endo | E | F | group_Endo |
| H1B_ES | G | H | group_ES |
>Tips
>1. A, B, C... mean the filenames of data, just like A.fastq.gz.
>2. If your data is .fastq.gz suffix file, please add the parameter of gzip, just like "--gzip true".
>3. If your filename of data is "Hela_cell_input.fastq.gz", please write its filename as "Hela_cell_input".

##### Comparefile
Comparefile is just like the following text which is a "\_vs\_" between two groups, just like the file [comparefile.txt](
https://github.com/kingzhuky/m6APipe/blob/master/test_data/comparefile.txt). 
>group_Endo_vs_group_ES


#### Pipeline Steps
m6APipe allows you to run pipelines skip the tools by your params.
You can skip the tools by using `--skip_ToolsName` or not(default), just list `--skip_metpeak`
| Step  | Pipeline |  Mode Parameter | Selection|
| :-: | :-: | :-: | :-: |
| Raw Data QC  | Fastp, FastQC                                                     |-|-|
| Reads Mapping      | star, bwa, tophat, hisat2                               |aligners|"star" OR "bwa" OR "tophat2" OR "hisat2" OR "none"|
| Sort BAM file AND Post-alignment QC  | samtools, RSeQC               |-|-|
| Peak Calling | MeTPeak, MACS2, MATK, meyer|peakCalling_mode |"group" OR "independence"|
| Combines Peaks information   | RobustRankAggreg, BEDtools, MSPC  |peakMerged_mode |"rank" OR "macs2" OR "MATK" OR "metpeak" OR "mspc"|
| Methylation analysis  | MeTDiff, QNB, MATK, Wilcox-test, DESeq2, edgeR |methylation_analysis_mode|"MATK" OR "QNB" OR "Wilcox-test" OR "MeTDiff" OR "edgeR" OR "DESeq2"|
| Expression analysis    | htseq, DESeq2, edgeR, cufflinks |expression_analysis_mode |"DESeq2" OR "edgeR" OR "none"|
### Dependencies

* Softwares
    * [Fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    * [Fastp](https://github.com/OpenGene/fastp)
    * [RSeQC](http://rseqc.sourceforge.net/)
    * [MultiQC](https://multiqc.info/)
    * [STAR](https://github.com/alexdobin/STAR)
    * [BWA](https://github.com/lh3/bwa)
    * [TopHat](https://ccb.jhu.edu/software/tophat/)
    * [HISAT2](https://ccb.jhu.edu/software/hisat2/)
    * [Bowtie2](https://github.com/BenLangmead/bowtie2)
    * [samtools](http://www.htslib.org/)
    * [MeTPeak](https://github.com/compgenomics/MeTPeak)
    * [MATK](http://matk.renlab.org)
    * [Meyer]()
    * [MACS2](https://github.com/taoliu/MACS)
    * [BEDtools](https://bedtools.readthedocs.io/en/latest/index.html)
    * [RobustRankAggreg](https://cran.r-project.org/web/packages/RobustRankAggreg/index.html)
    * [MeTDiff](https://github.com/compgenomics/MeTDiff)
    * [QNB](https://cran.r-project.org/src/contrib/Archive/QNB/)
    * [htseq](https://github.com/simon-anders/htseq)
    * [deseq2](http://bioconductor.org/packages/DESeq2/)
    * [edgeR](http://bioconductor.org/packages/edgeR/)
    * [cufflinks](http://cole-trapnell-lab.github.io/cufflinks/)
    * Several R packages for downstream analysis.

