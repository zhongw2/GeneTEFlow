################## BASE IMAGE ######################
FROM ubuntu:16.04

################## METADATA ######################
LABEL base_image="ubuntu:16.04"
LABEL about.summary="RNA-seq pipeline (Pfizer confidential)"
LABEL about.license="Pfizer"
LABEL about.tags="Genomics"


################## MAINTAINER ######################
MAINTAINER Xiaochuan (David) Liu <xiaochuan.liu@pfizer.com>

################## INSTALLATION ######################
#install tzdata for Time Zone Database
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y tzdata
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN dpkg-reconfigure -f noninteractive tzdata



#install the core libraries and software for ubunt
RUN  apt-get install --yes  default-jre perl python build-essential gcc-multilib   make gcc g++ zlib1g-dev zlib1g unzip zip wget


#install conda and update
RUN wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh -O Anaconda.sh
RUN bash Anaconda.sh -b -f -p /usr/local
RUN bash /usr/local/bin/activate
RUN rm Anaconda.sh
RUN conda update -n base conda
RUN conda config --system --add channels conda-forge
RUN conda config --system --append channels bioconda
RUN conda config --system --append channels r


#install software by conda
#Search tools listing on https://anaconda.org/
RUN conda install -c conda-forge r-base=3.5.1 
RUN conda install -c conda-forge r-devtools=1.13.6
RUN conda install -c conda-forge r-xlsx=0.6.1
RUN conda install -c conda-forge spyder=3.2.8
RUN conda install -c conda-forge r-lattice=0.20_35
RUN conda install -c conda-forge r-reshape2=1.4.3
RUN conda install -c conda-forge r-ggplot2=3.0.0
RUN conda install -c conda-forge r-rcolorbrewer=1.1_2
RUN conda install -c bioconda java-jdk=8.0.112

RUN apt-get update
#RUN conda update --all


#"LANG=C" is the default locale in the applications code itself,and get english-only in the menus.
#see https://ubuntuforums.org/showthread.php?t=1749224
ENV LANG=C


RUN mkdir -p  /RNASeq
COPY scripts  /RNASeq/scripts


