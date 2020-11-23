#!/bin/bash

mkdir -p ~/scratch/
echo export NXF_SINGULARITY_LOCALCACHEDIR="~/scratch/"
echo export NXF_SINGULARITY_CACHEDIR="~/scratch/"
echo export NXF_SINGULARITY_TMPDIR="~/scratch/"
echo export SINGULARITY_LOCALCACHEDIR="~/scratch/"
echo export SINGULARITY_CACHEDIR="~/scratch/"
echo export SINGULARITY_TMPDIR="~/scratch/"
echo export TMPDIR="~/scratch/"
echo export TEMPDIR="~/scratch/"
echo export SINGULARITY_DISABLE_CACHE=false

