FROM continuumio/miniconda3:24.1.2-0

LABEL Maintainer="Nicolai Meyerhoefer <nicolai.meyerhoefer@fau.de>"

USER root

RUN apt-get update

RUN conda install -y -c bioconda python=3.8 pandas=1.4.3 docopt=0.6.2 arboreto=0.1.6