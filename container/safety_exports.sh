#!/bin/bash

PREF=$2
USER=$1

sudo mkdir -p $PREF
sudo mkdir -p $PREF/nxf
sudo mkdir -p $PREF/sing
sudo mkdir -p $PREF/tmp
sudo chown -R $USER $PREF
sudo chgrp -R $USER $PREF

echo export NXF_SINGULARITY_LOCALCACHEDIR="$PREF/nxf"
echo export NXF_SINGULARITY_CACHEDIR="$PREF/nxf"
echo export NXF_SINGULARITY_TMPDIR="$PREF/nxf"

echo export SINGULARITY_LOCALCACHEDIR="$PREF/sing"
echo export SINGULARITY_CACHEDIR="$PREF/sing"
echo export SINGULARITY_TMPDIR="$PREF/sing"

echo export TMPDIR="$PREF/tmp"
echo export TEMPDIR="$PREF/tmp"

echo export SINGULARITY_ROOTFS="$PREF/sing"
echo export SINGULARITY_DISABLE_CACHE=false

