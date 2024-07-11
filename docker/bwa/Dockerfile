FROM ubuntu:24.04
ARG buildDate

LABEL about.image="Ubuntu base image to run bwa for sequence alignment"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.04"
LABEL image.version="0.1"
LABEL buildDate=${buildDate}
LABEL software="bwa to align fastq sequences"
LABEL software.version="0.7.18"
LABEL about.tags="Multiple Myeloma, ULP-WGS, bwa"

RUN apt-get update && apt-get install -y \
    bwa \
    zlib1g-dev 

ENTRYPOINT ["bwa"]