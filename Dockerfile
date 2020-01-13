FROM nfcore/base:1.7
LABEL authors="Kaiyu Zhu, Yu Sun" \
      description="Docker image containing all requirements for nf-core/meripseqpipe pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-meripseqpipe-1.0dev/bin:$PATH
