FROM r-base:4.2.2

LABEL Maintainer="Nicolai Meyerhoefer <nicolai.meyerhoefer@fau.de>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y gcc xml2 openssl libhdf5-dev libfftw3-dev curl

RUN R -e "install.packages(c('remotes', 'optparse', 'data.table', 'tidyverse', 'httr'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "remotes::install_github('satijalab/seurat', 'seurat5')"
