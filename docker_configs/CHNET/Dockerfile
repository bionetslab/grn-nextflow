FROM r-base:4.2.2

LABEL Maintainer="Nicolai Meyerhoefer <nicolai.meyerhoefer@fau.de>"

USER root

WORKDIR /

RUN apt-get update && apt-get install -y  gcc

RUN R -e "install.packages(c('remotes', 'optparse', 'data.table', 'igraph', 'stringr'), dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN R -e "remotes::install_github('Zhangxf-ccnu/chNet', subdir='pkg')"
