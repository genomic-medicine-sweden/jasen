#!/bin/bash

set -e

MNTROOT=/$(readlink -f .| cut -d"/" -f2)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ASSETS_DIR="${SCRIPT_DIR}/../.."
BIN_DIR="${ASSETS_DIR}/../bin"
TOKEN_DIR="${ASSETS_DIR}/.bigsdb_tokens"
CONTAINERS_DIR="${ASSETS_DIR}/../containers/"
CGMLST_DIR="${ASSETS_DIR}/cgmlst/klebsiella"

CLIENT_CREDENTIALS_FILE="${TOKEN_DIR}/client_credentials"

if [[ ! -f "$CLIENT_CREDENTIALS_FILE" ]]; then
    echo "Error: ${CLIENT_CREDENTIALS_FILE} not found. Copy assets/.bigsdb_tokens/client_credentials.template to jasen_dev/assets/.bigsdb_tokens/client_credentials and add your Pasteur API credentials." >&2
    exit 1
fi

if grep -q "insert_pasteur_client_id" "$CLIENT_CREDENTIALS_FILE"; then
    echo "Error: Replace the placeholder client_id value in ${CLIENT_CREDENTIALS_FILE} with your Pasteur API client_id." >&2
    exit 1
fi

if grep -q "insert_pasteur_client_secret" "$CLIENT_CREDENTIALS_FILE"; then
    echo "Error: Replace the placeholder client_secret value in ${CLIENT_CREDENTIALS_FILE} with your Pasteur API secret." >&2
    exit 1
fi

mkdir -p "${CGMLST_DIR}/alleles"

# Download alleles FASTA for all loci in scheme 18
apptainer exec --bind $MNTROOT ${CONTAINERS_DIR}/bactopia-py.sif \
    python ${ASSETS_DIR}/bin/bigsdb_downloader.py \
    --download_scheme \
    --key_name Pasteur \
    --site Pasteur \
    --token_dir "${TOKEN_DIR}" \
    --url "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18" \
    --output_dir "${CGMLST_DIR}/alleles"
