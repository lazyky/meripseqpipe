FROM nfcore/base
LABEL description="Docker image containing all requirements for nf-core/circpipe pipeline"

COPY environment.yml ./

ENV PATH /opt/conda/bin:$PATH
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/nf-core-m6APipe-1.0dev/bin:$PATH

#install mspc

