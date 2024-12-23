FROM rocker/verse:4.4.2
#Modified from here https://github.com/davismcc/r-tidybioc-img/blob/master/Dockerfile
LABEL description="Docker image containing DADA2, DECIPHER, phangorn, and latest R with tidyverse and Bioconductor packages"
 
# Install container-wide requrements gcc, pip, zlib, libssl, make, libncurses, fortran77, g++, R
RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get install -y --no-install-recommends \
        build-essential \
        git \
        libcurl4-gnutls-dev \
        libgsl-dev \
        libncurses5-dev \
        librdf0-dev \
        libreadline-dev \
        libssh2-1-dev \
        libssl-dev \
        libxml2-dev \
        libzmq3-dev \
        wget \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /usr/local/lib/R/site-library
ADD install.R /tmp/
RUN R -f /tmp/install.R

CMD ["R"]
