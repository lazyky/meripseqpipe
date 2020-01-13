From:nfcore/base:1.7
Bootstrap:docker

%labels
    DESCRIPTION Singularity image containing all requirements for the nf-core/meripseqpipe pipeline
    VERSION 1.0dev

%environment
    export PATH=/opt/conda/envs/nf-core-meripseqpipe-1.0dev/bin:$PATH
    export PATH=/opt/conda/bin:$PATH
    export PATH=/mspc:$PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
    /opt/conda/bin/conda install -y unzip 
    wget http://matk.renlab.org/download/MATK-1.0.jar
    wget https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz && \ 
    R CMD INSTALL QNB_1.1.11.tar.gz && \
    rm QNB_1.1.11.tar.gz
    git clone https://github.com/compgenomics/MeTDiff.git && \
    R CMD build MeTDiff/ && \
    R CMD INSTALL MeTDiff_1.0.tar.gz && \
    rm -rf MeTDiff*
    git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD build MeTPeak/ && \
    R CMD INSTALL MeTPeak_1.0.0.tar.gz && \
    rm -rf MeTPeak*
    wget -O mspc.zip "https://github.com/Genometric/MSPC/releases/download/v4.0.0/linux-x64.zip" && \
    unzip mspc.zip -d mspc && \
    chmod 775 mspc/mspc && \ 
    rm mspc.zip

