FROM rocker/r-ver:4.5
LABEL maintainer="Ryan Routsong <routsongrm@nih.gov>"
ENV DEBIAN_FRONTEND=noninteractive
ENV PACKAGES="libnss-sss gcc gfortran make libblas-dev liblapack-dev \
    libopenblas-dev liblapacke-dev wget ca-certificates libssl-dev libzstd-dev \
    libcurl4-openssl-dev libxml2-dev mono-mcs libncurses5-dev \
    libncursesw5-dev zlib1g zlib1g-dev git vim bzip2 libbz2-dev \
    liblzma-dev libhtscodecs2 locales curl locales-all libgmp3-dev"

RUN apt-get update && \
    apt-get install -y --no-install-recommends ${PACKAGES} && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8

ENV LANG=en_US.UTF-8 
ENV LC_ALL=en_US.UTF-8

# Write to the global Rprofile.site so settings persist in Singularity/Apptainer
# ($R_HOME/etc/Rprofile.site is inside the R install tree, not in $HOME)
RUN echo 'options(menu.graphics = FALSE)' >> "$(Rscript -e 'cat(R.home())')/etc/Rprofile.site" && \
    echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' >> "$(Rscript -e 'cat(R.home())')/etc/Rprofile.site"

COPY docker/install_packages.R /install_packages.R
RUN Rscript /install_packages.R && rm /install_packages.R

# make plink 1.9
RUN cd /opt; git clone https://github.com/chrchang/plink-ng.git
RUN cd /opt/plink-ng/1.9; ./plink_first_compile; make plink BLASFLAGS="-llapack -lopenblas"; ln -s /opt/plink-ng/1.9/plink /usr/bin/plink1.9
RUN cd /opt/plink-ng/2.0/build_dynamic; make; ln -s /opt/plink-ng/2.0/build_dynamic/plink2 /usr/bin/plink2.0
RUN mkdir -p /opt/plink-ng/2.3
COPY docker/plink2 /opt/plink-ng/2.3/plink2.3
RUN ln -s /opt/plink-ng/2.3/plink2.3 /usr/bin/plink2.3