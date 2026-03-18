FROM rocker/r-ver:4.5
LABEL maintainer="Ryan Routsong <routsongrm@nih.gov>"

ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="libnss-sss wget ca-certificates libssl-dev libcurl4-openssl-dev libxml2-dev \
    mono-mcs libncurses5-dev libncursesw5-dev zlib1g zlib1g-dev bzip2 libbz2-dev liblzma-dev libhtscodecs2 locales"

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

COPY docker/install_packages.R /install_packages.R
RUN Rscript /install_packages.R && rm /install_packages.R