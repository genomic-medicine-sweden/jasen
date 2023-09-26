#!/bin/bash

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
assdir="${scriptdir}/../"
containerdir="${scriptdir}/../../container/"
mntroot=/$(readlink -f .| cut -d"/" -f2)

cd $scriptdir
bash ./mlst-download_pub_mlst.sh
singularity exec --bind $mntroot ${containerdir}/blast.sif bash ${assdir}/mlst_db/mlst-make_blast_db.sh