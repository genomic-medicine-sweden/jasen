#!/bin/bash

scriptdir="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
assdir="${scriptdir}/../"
containers_dir="${scriptdir}/../../containers/"
mntroot=/$(readlink -f .| cut -d"/" -f2)

cd $scriptdir

# Download MLST data
if bash ./mlst-download_pub_mlst.sh; then
    echo "MLST download completed successfully"
else
    echo "MLST download failed" >&2
    exit 1
fi

# Create BLAST database
apptainer exec --bind $mntroot ${containers_dir}/blast.sif bash ${assdir}/mlst_db/mlst-make_blast_db.sh
