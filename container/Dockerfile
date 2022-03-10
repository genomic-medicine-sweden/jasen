FROM continuumio/miniconda3:latest AS condabuild
LABEL authors="Jyotirmoy Das" \
      description="Docker image containing all requirements for the JASEN pipeline"

COPY ./environment.yml /environment.yml
COPY ./environment.yaml /env.yml
RUN /opt/conda/bin/conda update conda && \
/opt/conda/bin/conda install mamba -c conda-forge && \
/opt/conda/bin/mamba env create -f /environment.yml && \
/opt/conda/bin/mamba env update -f /env.yml

FROM debian:buster-slim
RUN apt-get update && \
apt-get install -y git procps && \
apt-get clean -y
COPY --from=condabuild /opt/conda/envs/jasen /opt/conda/envs/jasen
ENV PATH=/opt/conda/envs/jasen/bin:$PATH
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
