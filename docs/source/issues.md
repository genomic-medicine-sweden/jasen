# Trouble shooting

Make sure that you have tried the following before submitting an issue.

* Always run the latest versions of the bioinformatical software.
* Verify you have execution permission for jasens `*.sif` images.
* Bind all the relevant paths (e.g. input, output, workdir) in docker/singularity/apptainer profile in nextflow.config, so they can be reached by containers 
* Old Singularity versions may sporadically produce the error `FATAL: could not open image jasen/container/*.sif: image format not recognised`

Please open an issue on the [issue tracker](https://github.com/genomic-medicine-sweden/jasen/issues) if you run into issues or bugs.
