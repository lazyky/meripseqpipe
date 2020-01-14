FROM nfcore/base:1.7
LABEL authors="Kaiyu Zhu, Yu Sun" \
      description="Docker image containing all requirements for nf-core/meripseqpipe pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-meripseqpipe-1.0dev/bin:$PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /mspc:$PATH

# install MATK
RUN wget http://matk.renlab.org/download/MATK-1.0.jar

# install QNB
RUN wget https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz && \ 
    R CMD INSTALL QNB_1.1.11.tar.gz && \
    rm QNB_1.1.11.tar.gz

# install MeTDiff
RUN git clone https://github.com/compgenomics/MeTDiff.git && \
    R CMD build MeTDiff/ && \
    R CMD INSTALL MeTDiff_1.0.tar.gz && \
    rm -rf MeTDiff*

# install MeTPeak
RUN git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD build MeTPeak/ && \
    R CMD INSTALL MeTPeak_1.0.0.tar.gz && \
    rm -rf MeTPeak*

# install MSPC
RUN conda install -y unzip 
RUN wget -O mspc.zip "https://github.com/Genometric/MSPC/releases/download/v4.0.0/linux-x64.zip" && \
    unzip mspc.zip -d mspc && \
    chmod 775 mspc/mspc && \ 
    rm mspc.zip
