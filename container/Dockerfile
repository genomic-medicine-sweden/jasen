FROM ubuntu:xenial
LABEL authors="Jyotirmoy Das" \
      description="Docker image containing all requirements for the JASEN pipeline"

RUN apt-get update && apt-get install -y \
 wget \
 git \
 --no-install-recommends -y build-essential make curl unzip locales \
 build-essential && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/* 

RUN wget --no-check-certificate https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
 bash ./Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b  && \
 rm Miniconda3-latest-Linux-x86_64.sh

ENV PATH="/miniconda/bin:$PATH"
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

COPY container/main.yml /environment.yml
COPY container/reporting.yml /env.yml
RUN conda update conda && \
conda env create -f /environment.yml && \
conda env update -f /env.yml && \
echo "source activate jasen" > /etc/bash.bashrc
ENV PATH="/miniconda/envs/jasen/bin:$PATH"

CMD . activate jasen
