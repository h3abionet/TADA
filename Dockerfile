FROM bioconductor/release_core2

LABEL project="dada2"
LABEL version="release"

MAINTAINER Katie Lennard

COPY --from: dada2_procps:latest

RUN R -e 'source("https://bioconductor.org/biocLite.R"); biocLite("dada2"); biocLite("DECIPHER"); biocLite("biomformat")'
RUN R -e 'install.packages(c("phangorn","dplyr"), dependencies=TRUE)'

CMD ["R"]

###########################################################
RUN apt-get update && apt-get install -y procps

################## Hex specific ###########################
RUN mkdir -p /researchdata/fhgfs
