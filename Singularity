From:nfcore/base
Bootstrap:docker

%labels
    DESCRIPTION Singularity image containing all requirements for the nf-core/m6APipe pipeline
    VERSION 1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-m6APipe-1.0dev/bin:$PATH
    export PATH
    export PATH=$PATH:/home/wqj/miniconda3/bin
    export PATH=$PATH:/home/wqj/miniconda3/envs/tools_in_python3/bin
    export PATH=$PATH:/home/wqj/miniconda3/envs/tools_in_python2/bin

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
    conda install --yes nextflow fastqc bedtools ucsc-gtftogenepred hiast2 bowtie2 bwa samtools star tophat
    conda install --yes bioconductor-edger bioconductor-deseq2 htseq cufflinks bioconductor-exomepeak macs2 python=2.7.13
    conda install --yes -c bioconda rseqc


