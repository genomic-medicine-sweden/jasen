#!/usr/bin/env bash

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

NAME=$1
if [ $NAME ]; then
    conda env create --name $NAME -f $scriptdir/environment.yaml
else
    conda env create -f $scriptdir/environment.yaml
fi
