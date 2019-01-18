# nf-core/m6APipe: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the RSeQC report, which summarises results at the end of the pipeline.


## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [RSeQC](#RSeQC) - evaluate high throughput sequence data especially RNA-seq data

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the RSeQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images


## RSeQC
[RSeQC](http://rseqc.sourceforge.net/) provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data. Some basic modules quickly inspect sequence quality, nucleotide composition bias, PCR bias and GC bias, while RNA-seq specific modules evaluate sequencing saturation, mapped reads distribution, coverage uniformity, strand specificity, transcript level RNA integrity etc.

The pipeline has special steps which allow the software versions used to be reported in the RSeQC output for future traceability.

**Output directory: `results/RSeQC`**

* `Project_RSeQC_report.html`
  * RSeQC **
* `Project_RSeQC_data/`
  * Direct****

For more information about how to use RSeQC reports, see http://rseqc.sourceforge.net/
