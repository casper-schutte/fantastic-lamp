FROM continuumio/miniconda3

ADD environment.yaml environment.yaml
RUN conda env update -n base --file environment.yaml --prune && conda clean -faipy
ENV LD_PRELOAD=/opt/conda/lib/libjemalloc.so.2
ADD . /repo/
