FROM continuumio/anaconda3

RUN apt-get update
RUN apt-get -y install r-base-core

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

WORKDIR /opt
RUN mkdir -p minikraken
WORKDIR /opt/minikraken
# UPDATE KRAKEN DB FOR NEWER VERSIONS
# MOVING AWAY FROM KRAKEN2 - POORER FUNCTIONALITY
# RUN wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken2_v2_8GB_201904_UPDATE.tgz
RUN wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171101_4GB_dustmasked.tgz
RUN tar -xzvf minikraken_20171101_4GB_dustmasked.tgz
RUN mv minikraken_20171101_4GB_dustmasked minikraken

## SWITCH TO MAMBA
RUN conda install -c conda-forge mamba

RUN mamba config --add channels conda-forge
RUN mamba config --add channels bioconda
RUN mamba install -y kraken perl biopython shovill coverm

# Install R packages
RUN Rscript -e "install.packages('ggplot2',lib='/usr/lib/R/library',dependencies=TRUE,repos='http://cran.rstudio.com')"

# Set up environment for coverage estimation
RUN mamba create -n coverage coverm

# Set up environment for QC
RUN mamba create -n qc fastqc multiqc

# REMEMBER TO PULL LATEST CHANGES FROM GIT REPOSITORY
RUN cp -r *.py /usr/bin
RUN cp -r *.sh /usr/bin
RUN cp -r *.R /usr/bin

WORKDIR /mnt
RUN useradd -ms /bin/bash -g sudo scrub
USER scrub

LABEL version="1.4"
MAINTAINER Ola Brynildsrud "olbb@fhi.no" 
