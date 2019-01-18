# Example usage for non-human species 

## Introduction 
LncPipe accepts raw reads, annotations and genome reference as input to conduct the whole analysis, 
and is also applicable for selected referenced species. In the first version, we mainly focus on human 
because of well-organized lncRNA annotation files (with .gtf suffixed) from GENCODE and LNCipedia. 
For non-human species, user are required to provide both protein coding annotation file and lncRNA annotation file separately, 
which could be download from GENCODE(human or mouse only) or Ensemble databases. However, not all non-human species are supported 
at present, since one essential tool CPAT included in LncPipe only available for 4 species (human, mouse, fly and zebrafish). 

## Example usage for mouse 
### Step 1. Prepare input files. The following files required by LncPipe
* hisat index: ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz
```shell
    #if database not exsit 
    mkdir ~/database/mouse
    cd ~/database/mouse  
    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grcm38_tran.tar.gz
    tar -xzvf grcm38_tran.tar.gz
```
* Genome reference: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz
```shell
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/GRCm38.p5.genome.fa.gz
    gunzip GRCm38.p5.genome.fa.gz
```

* GTF files including all features :ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
```shell
    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M16/gencode.vM16.annotation.gtf.gz
    gunzip gencode.vM16.annotation.gtf.gz
```

* Raw sequence file with \*.fastq.gz / \*.fastq suffixed
> For example there are 4 samples with eight .gz suffixed fastq files, like
```shell
sample1_control_input_1.fastq.gz,
sample1_control_ip_1.fastq.gz,
sample1_treated_input_1.fastq.gz,
sample1_treated_ip_1.fastq.gz,
sample2_control_input_1.fastq.gz,
sample2_control_ip_1.fastq.gz,
sample2_treated_input_1.fastq.gz,
sample2_treated_ip_1.fastq.gz,
```

* `Design.file` are required if you are going to perform differential expression analysis between groups. 

### Step 2. Edit `nextflow.config` file 
Leave the other line unchanged, modified the following sentences like below (According to the location of download files):
> If your are using docker in your server, plz modify `docker.config` instead.  

```shell
    fastq_ext = '*_{1,2}.fastq.gz'
    fasta_ref = ''
    design = 'designfile.txt'
    hisat2_index = ''
    mspcpath='/opt/mspc_v3.3'
    //human gtf only
    gencode_annotation_gtf = ""
    lncipedia_gtf = null
    species="mouse"// mouse , zebrafish, fly
    known_coding_gtf="gencode.vM16.proteincoding.gtf"
    known_lncRNA_gtf="gencode.vM16.annotation.gtf"

```
### Step 3. Start your analysis trip with command below   

```shell
    nextflow run -with-trace -with-report report.html -with-timeline timeline.html main.nf 
    #or running in a docker image  
    nextflow -c docker.config  main.nf 
```
> The default running tools in each step are fastp, hisat, gffcompare, stringtie, cpat, plek, sambamba, kallisto ,edgeR and LncPipeReporter, if you want to change the tool in each step, plz modify `config` file instead.

* Any question, plz open an issue in the issue page, we will reply ASAP :)
