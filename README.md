# m6APipe
**MeRIP-seq analysis pipeline**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker / singularity containers making installation trivial and results highly reproducible.


### Documentation   

A full tutorial of m6APipe can be found at Wiki page of this project. plz go to the https://github.com/kingzhuky/m6APipe/wiki


### Pipeline Steps

m6APipe allows you to run pipelines skip the tools by your params.
You can skip the tools by using `--skip_ToolsName` or not(default).


| Step                                    | Pipeline                        | Pipeline(skip_mode)             |
|-----------------------------------------|---------------------------------|---------------------------------|
| Raw Data QC                             | FastQC                          |----  --skip_aligners true  -----|
| Reads Mapping                           | star, bwa, tophat, hisat2       |----  --skip_aligners true  -----|
| Sort BAM file AND Post-alignment QC     | samtools, RSeQC                 |---------------------------------|
| Reads counting                          | htseq-count                     |---------------------------------|
| Peak Calling                            | MeTPeak, exomePeak, macs2, MATK |---  --skip_peakCalling true  ---|
| Differential methylation analysis       | MeTDiff, exomePeak, QNB, MATK   |-  --skip_diffpeakCalling true  -|
| Differential expression analysis        | deseq2, edgeR, cufflinks        |---  --skip_expressiion true  ---|
| Combines Peaks information              | MSPC                            |------  --skip_mspc true  -------|

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
    * [exomePeak](http://bioconductor.org/packages/exomePeak/)
    * [macs](https://github.com/taoliu/MACS)
    * [MeTDiff](https://github.com/compgenomics/MeTDiff)
    * [QNB](https://cran.r-project.org/src/contrib/Archive/QNB/)
    * [MSPC](https://github.com/Genometric/MSPC)

    * Several R packages for downstream analysis.