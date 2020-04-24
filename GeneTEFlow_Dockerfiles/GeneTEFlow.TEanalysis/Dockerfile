################## BASE IMAGE ######################
#FROM ubuntu:16.04
FROM	continuumio/miniconda3
################## METADATA ######################
#LABEL base_image="ubuntu:16.04"
LABEL base_image="/miniconda3"
LABEL about.summary="locus-specific regulation of interspersed repeat expression"
#LABEL about.summary="RNA-seq pipeline (Pfizer confidential)"
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


RUN    conda install  --override-channels -c iuc -c bioconda -c conda-forge -c defaults -c r python=2.7.13 bioconductor-deseq2=1.16.1 r-base=3.4.1 r-pheatmap bioconductor-vsn bioconductor-biocparallel=1.12.0 r-ggrepel star=2.5.3a bedtools=2.25.0 samtools=1.1 stringtie=1.3.3 igvtools=2.3.93 ucsc-genepredtobed ucsc-gtftogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig r-hexbin

#RUN	conda create --name squire --override-channels -c iuc -c bioconda -c conda-forge -c defaults -c r python=2.7.13 bioconductor-deseq2=1.16.1 r-base=3.4.1 r-pheatmap bioconductor-vsn bioconductor-biocparallel=1.12.0 r-ggrepel star=2.5.3a bedtools=2.25.0 samtools=1.1 stringtie=1.3.3 igvtools=2.3.93 ucsc-genepredtobed ucsc-gtftogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig r-hexbin 
#RUN	echo "source activate squire" > ~/.bashrc

#ENV PATH /opt/conda/envs/squire/bin:$PATH

RUN	git clone https://github.com/wyang17/SQuIRE 

WORKDIR /SQuIRE
RUN	cd /SQuIRE 
RUN	pip install -e .
#RUN	cd SQuIRE \
#RUN	python setup.py install
RUN	cd /SQuIRE 
RUN	git pull
#ENV PATH /SQuIRE/squire:$PATH
COPY  Count.md.py  /SQuIRE/squire/Count.py

RUN apt-get install -y  --reinstall procps 





#RUN apt-get update
#RUN conda update --all


#"LANG=C" is the default locale in the applications code itself,and get english-only in the menus.
#see https://ubuntuforums.org/showthread.php?t=1749224




