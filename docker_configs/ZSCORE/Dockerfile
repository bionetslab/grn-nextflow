FROM r-base:4.2.2

LABEL Maintainer="Nicolai Meyerhoefer <nicolai.meyerhoefer@fau.de>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y  gcc

RUN R -e "install.packages(c('BiocManager', 'optparse', 'data.table'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "BiocManager::install('dcanr')"
