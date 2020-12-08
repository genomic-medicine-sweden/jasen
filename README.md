# JASEN
_Json producing Assembly driven microbial Sequence analysis pipeline to support Epitypification and Normalize classification decisions_

## Setup
* git clone --recurse-submodules --single-branch --branch master  https://github.com/genomic-medicine-sweden/JASEN.git
* Edit `nextflow.config`
* _`Optionally: bash safety_exports.sh USER PREFIX`_


### Singularity usage
* Install Singularity
* `bash container/build_container.sh`
* `singularity exec -B JASEN_INSTALL_DIR:/external -B WORKDIR:/out IMAGE nextflow -C /external/nextflow.config run /JASEN/main.nf -profile local,singularity`


### Conda usage
* Install Conda ( https://www.anaconda.com/distribution )
* Install nextFlow ( `curl -s https://get.nextflow.io | bash` )
* `bash setup.sh`
* `nextflow run main.nf -profile -local,conda`
