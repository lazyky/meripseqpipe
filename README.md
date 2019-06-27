# m6APipe
**MeRIP-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible. N6-methyladenosine (m6A) is the most prevalent modification in the mRNA of many eukaryotic species, including yeast, plants, flies, and mammals. In order to analyze m6A-seq data, we developed a user-friendly, integrated analysis pipeline called m6APipe based on Nextflow. It integrated ten main functional modules including preprocessing, QC, read mapping, peak calling, merging peaks, differential methylation analysis, differential expression analysis, motif search, annotation, and data visualization. 


### Documentation   

A full tutorial of m6APipe can be found at Wiki page of this project. plz go to the https://github.com/kingzhuky/m6APipe/wiki


### Pipeline Steps

m6APipe allows you to run pipelines skip the tools by your params.
You can skip the tools by using `--skip_ToolsName` or not(default).

| Step                                    | Pipeline                        | Pipeline(skip_mode)             |
|-----------------------------------------|---------------------------------|---------------------------------|
| Raw Data QC                             | Fastp, FastQC                   |---------------------------------|
| Reads Mapping                           | star, bwa, tophat, hisat2       |------  --aligners "star"  ------|
| Sort BAM file AND Post-alignment QC     | samtools, RSeQC                 |---------------------------------|
| Peak Calling                            | MeTPeak, MACS2, MATK, meyer     |---  --skip_peakCalling true  ---|
| Combines Peaks information              | RobustRankAggreg, BEDtools      |---------------------------------|
| Methylation analysis                    | MeTDiff, QNB, MATK, Wilcox-test |-  --skip_diffpeakCalling true  -|
| Expression analysis                     | htseq, deseq2, edgeR, cufflinks |---  --skip_expressiion true  ---|


### Dependencies
* Softwares
    * [Fastqc](https://github.com/OpenGene/fastp)
    * [STAR](https://github.com/alexdobin/STAR)
    * [BWA](https://github.com/lh3/bwa)
    * [TopHat](https://ccb.jhu.edu/software/tophat/)
    * [HISAT2](https://ccb.jhu.edu/software/hisat2/)
    * [Bowtie2](https://github.com/BenLangmead/bowtie2)
    * [samtools](http://www.htslib.org/)
    * [htseq](https://github.com/simon-anders/htseq)
    * [deseq2](http://bioconductor.org/packages/DESeq2/)
    * [edgeR](http://bioconductor.org/packages/edgeR/)
    * [MeTPeak](https://github.com/compgenomics/MeTPeak)
    * [macs](https://github.com/taoliu/MACS)
    * [MeTDiff](https://github.com/compgenomics/MeTDiff)
    * [QNB](https://cran.r-project.org/src/contrib/Archive/QNB/)

    * Several R packages for downstream analysis.