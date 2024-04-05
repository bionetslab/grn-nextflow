FROM continuumio/miniconda3:24.1.2-0

LABEL Maintainer="Nicolai Meyerhoefer <nicolai.meyerhoefer@fau.de>"

USER root

RUN apt-get update && apt-get install -y gcc

RUN conda install -y -c bioconda python=3.8 pandas numpy docopt matplotlib cython pydot graphviz scipy 

RUN git clone https://github.com/gihannagalindez/boostdiff_inference.git

WORKDIR /boostdiff_inference/

RUN pip install .

WORKDIR /
