#!/usr/bin/env bash

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

NAME=$1
if [ $NAME ]; then
    conda env create --name $NAME -f $scriptdir/requirements.txt
    $scriptdir/../..//bin/pipeline_result_processor
    pip install .
    cd $scriptdir
else
    conda env create --name jasen -f $scriptdir/requirements.txt
    $scriptdir/../..//bin/pipeline_result_processor
    pip install .
    cd $scriptdir
fi
