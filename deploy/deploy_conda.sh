#!/usr/bin/env bash

NAME=$1
scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

if ! [ $NAME ]; then
   NAME=jasen
fi   


#Unload environment
conda info | tac | tac | grep -q $NAME && source deactivate || :
#Remove environment if already present
echo "Purging any duplicate existing environment"
conda remove -y -n $NAME --all || :

echo "Creating JASEN environment named $NAME"
conda env create --name $NAME -f $scriptdir/reqs/requirements.txt -q 
source activate jasen
cd $scriptdir/../bin/pipeline_result_processor
pip install .
cd $scriptdir
echo "Done!"
