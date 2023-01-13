#!/usr/bin/env bash

# quit script if any of the commands fails. Note that && commands should be in paranthesis for this to work.
set -eox pipefail
trap 'exit_status="$?" && echo Failed on line: $LINENO at command: $BASH_COMMAND && echo "exit status $exit_status" && exit' ERR

# name for the conda environment. Should be possible to create another name instead of jasen.
NAME=$1
# script directory
scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# If the user did not gave a name for the conda environment on the command line, then the name will be jasen
if ! [ $NAME ]; then
   NAME=jasen
fi

conda info
#Unload environment. Deactivation not working properly.
conda info | grep -q $NAME && conda deactivate || :
conda info

#Remove environment if already present
echo "Purging any duplicate existing environment"
conda remove -y -n $NAME --all || :

echo "Creating JASEN environment named $NAME"
#conda create --name $NAME -f $scriptdir/reqs/reqs.txt -q 
# mamba is used for quicker env resolve time
mamba env create -f deploy/reqs/env-index-mini.yaml 
source activate jasen
conda info
cd $scriptdir/../bin/pipeline_result_processor
pip install .
cd $scriptdir
echo "Done!"
