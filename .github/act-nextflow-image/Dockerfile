# Docker image to be used for running CI locally using ACT
FROM ghcr.io/catthehacker/ubuntu:act-latest

# Install Java and Singularity
RUN apt-get update \
    && sudo apt-get install -y openjdk-21-jre \
        libfuse2 \
        crun \
        runc \
        cryptsetup-bin \
        fuse \
        fuse2fs \
        squashfs-tools \
        uidmap \
    && wget https://github.com/sylabs/singularity/releases/download/v4.1.0/singularity-ce_4.1.0-jammy_amd64.deb \
    && sudo dpkg -i singularity-ce_4.1.0-jammy_amd64.deb

# Install Nextflow
RUN curl -fsSL get.nextflow.io | bash && sudo mv nextflow /usr/local/bin
