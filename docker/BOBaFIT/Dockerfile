FROM r-base:4.4.1
ARG buildDate

ADD /BOBaFIT-lib/ /BOBaFIT/

LABEL about.image="R base image to run BOBaFIT for diploid region correction on seg files and chromosome arm level categorization"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="r-base:4.4.1"
LABEL image.version="0.1"
LABEL buildDate=${buildDate}
LABEL software="BOBaFIT to diploid region correction of seg files"
LABEL software.version="1.8.0"
LABEL about.tags="Multiple Myeloma, ULP-WGS, BOBaFIT"

RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libssh2-1-dev \
    procps 
RUN R -e 'install.packages(c("BiocManager", "data.table", "dplyr", "magrittr", "NbClust","optparse", "stringr"))'
RUN R -e 'BiocManager::install("GenomicRanges")'

ENTRYPOINT ["Rscript"]