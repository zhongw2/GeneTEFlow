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
RUN conda config --system --append channels r


#install software by conda
#Search tools listing on https://anaconda.org/
RUN conda install -c conda-forge r-base=3.4.1
RUN conda install -c conda-forge r-devtools=1.13.6
RUN conda install -c conda-forge spyder=3.2.8
RUN conda install -c conda-forge r-reshape2=1.4.3
RUN conda install -c conda-forge r-ggplot2=3.0.0
RUN conda install -c conda-forge r-gridextra=2.3
RUN conda install -c conda-forge r-gridbase=0.4_7
RUN conda install -c conda-forge r-gtable=0.2.0
RUN conda install -c conda-forge r-rcolorbrewer=1.1_2
RUN conda install -c conda-forge r-gplots=3.0.1
RUN conda install -c bioconda bioconductor-complexheatmap=1.18.1
RUN conda install -c conda-forge r-circlize=0.4.4
RUN conda install -c conda-forge r-venndiagram=1.6.20
RUN conda install -c conda-forge r-survival=2.42_6
RUN conda install -c conda-forge r-formula=1.2_3
RUN conda install -c conda-forge r-latticeextra=0.6_28
RUN conda install -c conda-forge r-cluster=2.0.7_1
RUN conda install -c conda-forge r-rpart=4.1_13
RUN conda install -c conda-forge r-nnet=7.3_12
RUN conda install -c conda-forge r-acepack=1.4.1
RUN conda install -c conda-forge r-foreign=0.8_71
RUN conda install -c conda-forge r-data.table=1.11.4
RUN conda install -c conda-forge r-htmltable=1.12
RUN conda install -c conda-forge r-viridislite=0.3.0
RUN conda install -c conda-forge r-htmltools=0.3.6
RUN conda install -c conda-forge r-base64enc=0.1_3
RUN conda install -c conda-forge r-hmisc=4.1_1
RUN conda install -c conda-forge r-knitr=1.20

#RUN apt-get update
#RUN conda update --all


#install tables_0.8.4 by R
RUN wget https://cran.r-project.org/src/contrib/tables_0.8.8.tar.gz
RUN R CMD INSTALL tables_0.8.8.tar.gz
RUN rm tables_0.8.8.tar.gz



#"LANG=C" is the default locale in the applications code itself,and get english-only in the menus.
#see https://ubuntuforums.org/showthread.php?t=1749224
ENV LANG=C


RUN mkdir -p  /RNASeq
COPY ngsdb  /RNASeq/ngsdb
COPY scripts  /RNASeq/scripts
COPY apps  /RNASeq/apps




