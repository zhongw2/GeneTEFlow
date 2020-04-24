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

#install the core libraries and software for ubuntu
RUN  apt-get install --yes  default-jre perl python build-essential gcc-multilib   make gcc g++ zlib1g-dev zlib1g unzip zip wget

#install conda and update
RUN wget https://repo.anaconda.com/archive/Anaconda2-5.2.0-Linux-x86_64.sh -O Anaconda.sh
RUN bash Anaconda.sh -b -f -p /usr/local
RUN bash /usr/local/bin/activate
RUN rm Anaconda.sh
RUN conda update -n base conda
RUN conda config --system --add channels conda-forge
RUN conda config --system --append channels bioconda


#install software by conda
#Search tools listing on https://anaconda.org/
RUN conda install -c conda-forge r-base=3.4.1
RUN conda install -c conda-forge r-devtools=1.13.5
RUN conda install -c bioconda fastqc=0.11.7
RUN conda install -c bioconda hisat2=2.1.0
RUN conda install -c bioconda samtools=1.3.1
RUN conda install -c bioconda stringtie=1.3.4
RUN conda install -c bioconda star=2.6.0c
RUN conda install -c bioconda rsem=1.3.0
RUN conda install -c bioconda bioconductor-deseq2=1.18.1
RUN conda install -c conda-forge libiconv=1.15
RUN conda install -c bioconda bioconductor-tximport=1.8.0



RUN conda install -c conda-forge r-readr=1.1.1
RUN conda install -c bioconda bioconductor-genomeinfodbdata=1.1.0 
  
RUN conda install -c conda-forge r-plyr=1.8.4
RUN conda install -c bioconda trimmomatic=0.36

#"LANG=C" is the default locale in the applications code itself,and get english-only in the menus.
#see https://ubuntuforums.org/showthread.php?t=1749224
ENV LANG=C

# Set the locale
RUN apt-get install -y locales
RUN locale-gen en_US.UTF-8 

RUN mkdir -p  /RNASeq
COPY ngsdb  /RNASeq/ngsdb
COPY scripts  /RNASeq/scripts


