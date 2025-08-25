#!/bin/bash

set -e

MNTROOT=/$(readlink -f .| cut -d"/" -f2)
ASSETS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
CONTAINERS_DIR="${ASSETS_DIR}/../containers/"

CLIENT_ID=""
CLIENT_SECRET=""

PUBMLST_ID=${CLIENT_ID:=$PUBMLST_CLIENT_ID}
PUBMLST_SECRET=${CLIENT_SECRET:=$PUBMLST_CLIENT_SECRET}

# Create PubMLST access token
apptainer exec --bind $MNTROOT ${CONTAINERS_DIR}/bactopia-py.sif bactopia-pubmlst-setup \
    --client-id $PUBMLST_ID \
    --client-secret $PUBMLST_SECRET \
    -d pubmlst_saureus_seqdef \
    -sd $ASSETS_DIR

# Create PubMLST database
apptainer exec --bind $MNTROOT ${CONTAINERS_DIR}/bactopia-py.sif bactopia-pubmlst-build \
    --force \
    -d all \
    -t $ASSETS_DIR \
    -o $ASSETS_DIR
