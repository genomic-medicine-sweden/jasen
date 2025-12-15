#!/bin/bash

set -e

MNTROOT=/$(readlink -f .| cut -d"/" -f2)
ASSETS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
CONTAINERS_DIR="${ASSETS_DIR}/../containers/"

CLIENT_ID=""
CLIENT_SECRET=""

PUBMLST_ID=${CLIENT_ID:=$PUBMLST_CLIENT_ID}
PUBMLST_SECRET=${CLIENT_SECRET:=$PUBMLST_CLIENT_SECRET}

# Check if required environment variables are set
if [[ -z "$PUBMLST_ID" ]]; then
    echo "Error: PUBMLST_CLIENT_ID environment variable is not set or CLIENT_ID is not defined in this script" >&2
    exit 1
fi

if [[ -z "$PUBMLST_SECRET" ]]; then
    echo "Error: PUBMLST_CLIENT_SECRET environment variable is not set or CLIENT_SECRET is not defined in this script" >&2
    exit 1
fi

# Create PubMLST access token
apptainer exec --bind $MNTROOT ${CONTAINERS_DIR}/bactopia-py.sif bactopia-pubmlst-setup \
    --force \
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
