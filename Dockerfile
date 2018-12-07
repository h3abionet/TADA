FROM bioconductor/release_core2

LABEL project="dada2"
LABEL version="release"

MAINTAINER Katie Lennard

RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("dada2"); biocLite("DECIPHER"); biocLite("biomformat")'
RUN R -e 'install.packages(c("phangorn","dplyr"), dependencies=TRUE)'

CMD ["R"]

###########################################################
FROM debian:jessie
ENV CONDA_INSTALLER="Miniconda3-latest-Linux-x86_64.sh"
#Exports conda path
ENV PATH $PATH:/opt/conda/bin/
RUN apt-get update && apt-get install -y procps

################## Hex specific ###########################
RUN mkdir -p /researchdata/fhgfs
