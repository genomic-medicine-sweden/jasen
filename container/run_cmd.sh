
singularity exec -B /home/is/work/code/JASEN/:/external -B /tmp/:/out jasen_2020-12-02.sif nextflow -C /external/nextflow.config run /JASEN/main.nf -profile local,singularity
