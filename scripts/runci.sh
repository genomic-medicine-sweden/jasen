#!/bin/bash
# This script uses the act tool to execute GitHub Actions workflows locally
# for increased iteration speed.
# It requires that you have Docker installed.
# It requires the act [1] binary to be available in your path.
# It also requires a late version of the GitHub commandline client (`gh`) [2]
# to be installed and available on your path.
# [1] https://github.com/nektos/act
# [2] https://github.com/cli/cli

# This enables you to send extra options to act, such as
# --pull=false to avoid pullimg images from Dockerhub.
# or the -v option to get more output for troubleshooting.

act_options=$1;
act \
    -s GITHUB_TOKEN="$(gh auth token)" \
    $act_options \
    -e scripts/event.json \
    -P ubuntu-latest=genomicmedicinesweden/act-nextflow:latest \
    -W .github/workflows/ci.yml \
|& tee ci-output-$(date +%Y%m%d.%H%M%S).log
