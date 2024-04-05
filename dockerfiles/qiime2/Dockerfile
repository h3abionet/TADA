FROM qiime2/core
MAINTAINER Gerrit Botha "gerrit.botha@uct.ac.za"

RUN apt-get -qq update && apt-get -qq -y install procps
ENV TZ=Africa/Johannesburg
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
