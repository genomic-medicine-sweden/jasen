#!/usr/bin/env bash

NAME=$1
scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if ! [ $NAME ]; then
   NAME=jasen
fi   


#Unload environment
conda info | grep -q $NAME && source deactivate || :
#Remove environment if already present
echo "Purging any duplicate existing environment"
conda remove -y -n $NAME --all || :

echo "Creating JASEN environment named $NAME"
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -y -n $NAME mamba && conda activate $NAME
mamba install -y ncbi-amrfinderplus bwa blast chewbbaca
#mamba env create -f deploy/reqs/env-index-mini.yaml 
echo "Done!"
