FROM rocker/r-ver:4.5
LABEL maintainer="Ryan Routsong <routsongrm@nih.gov>"

ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="libnss-sss wget ca-certificates libssl-dev libcurl4-openssl-dev libxml2-dev \
    mono-mcs libncurses5-dev libncursesw5-dev zlib1g zlib1g-dev bzip2 libbz2-dev liblzma-dev libhtscodecs2 locales"

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Write to the global Rprofile.site so settings persist in Singularity/Apptainer
# ($R_HOME/etc/Rprofile.site is inside the R install tree, not in $HOME)
RUN echo 'options(menu.graphics = FALSE)' >> "$(Rscript -e 'cat(R.home())')/etc/Rprofile.site" && \
    echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> "$(Rscript -e 'cat(R.home())')/etc/Rprofile.site"

COPY docker/install_packages.R /install_packages.R
RUN Rscript /install_packages.R && rm /install_packages.R