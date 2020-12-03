#!/bin/bash



sudo mkdir -p /scratch/
sudo chown is /scratch/
mkdir -p /scratch/nxf
mkdir -p /scratch/sing
mkdir -p /scratch/tmp

echo export NXF_SINGULARITY_LOCALCACHEDIR="/scratch/nxf"
echo export NXF_SINGULARITY_CACHEDIR="/scratch/nxf"
echo export NXF_SINGULARITY_TMPDIR="/scratch/nxf"

echo export SINGULARITY_LOCALCACHEDIR="/scratch/sing"
echo export SINGULARITY_CACHEDIR="/scratch/sing"
echo export SINGULARITY_TMPDIR="/scratch/sing"

echo export TMPDIR="/scratch/tmp"
echo export TEMPDIR="/scratch/tmp"

echo export SINGULARITY_ROOTFS="/scratch/sing"
echo export SINGULARITY_DISABLE_CACHE=false

