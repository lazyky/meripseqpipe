FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/m6APipe pipeline"

COPY environment.yml ./

ENV PATH /opt/conda/bin:$PATH
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-m6APipe-1.0dev/bin:$PATH

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