#!/bin/bash

set -e

MNTROOT=/$(readlink -f .| cut -d"/" -f2)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ASSETS_DIR="${SCRIPT_DIR}/../.."
TOKEN_DIR="${ASSETS_DIR}/.bigsdb_tokens"
CONTAINERS_DIR="${ASSETS_DIR}/../containers/"
CGMLST_DIR="${ASSETS_DIR}/cgmlst/klebsiella"

ACCESS_TOKENS_FILE="${TOKEN_DIR}/access_tokens"

if [[ ! -f "$ACCESS_TOKENS_FILE" ]]; then
    echo "Error: ${ACCESS_TOKENS_FILE} not found. Create it and add your Pasteur API credentials." >&2
    exit 1
fi

if grep -q "insert_pasteur_token_here" "$ACCESS_TOKENS_FILE"; then
    echo "Error: Replace the placeholder token value in ${ACCESS_TOKENS_FILE} with your Pasteur API token." >&2
    exit 1
fi

if grep -q "insert_pasteur_secret_here" "$ACCESS_TOKENS_FILE"; then
    echo "Error: Replace the placeholder secret value in ${ACCESS_TOKENS_FILE} with your Pasteur API secret." >&2
    exit 1
fi

mkdir -p "${CGMLST_DIR}/alleles"

# Download Klebsiella cgMLST schema from BIGSdb Pasteur
apptainer exec --bind $MNTROOT ${CONTAINERS_DIR}/bactopia-py.sif \
   bigsdb_downloader.py \
   --key_name Pasteur \
   --site Pasteur \
   --token_dir "${TOKEN_DIR}" \
   --url "https://bigsdb.pasteur.fr/api/db/pubmlst_klebsiella_seqdef/schemes/18/alleles_fasta" \
   --output_dir "${CGMLST_DIR}/alleles" \
   --cron
