# Installation

## Requirements

* Singularity
* Nextflow (`curl -s https://get.nextflow.io | bash`)

**Recommended**

* Conda
* Singularity Remote Login

## Development deployment (self-contained)

### Copy code locally

```bash
git clone --recurse-submodules --single-branch --branch master \\ 
    git@github.com:genomic-medicine-sweden/jasen.git &&        \\
cd jasen
```

### Access to OCI registries (Optional)

```bash
singularity remote login
```

### Create singularity images. 

The containers will be attempted to be built and downloaded as part of the main Makefile (that is, when running `make install` in the main repo folder).

```bash
cd container
make
```

### Download references and databases using singularity. 

First, make sure you stand in the main jasen folder (so if you cd:ed into the `container` folder before, you need to cd back to the main folder with `cd ..`). Then run the `install` make rule:

```bash
make install
```

Finally, run checks:

```bash
make check
```

Any errors produced during this step will hinder pipeline execution in unexpected ways.

## Configuration and test data

### Config 

Source: `configs/nextflow.base.config`

* Edit the `root` parameter in `configs/nextflow.base.config`
* Edit the `krakenDb`, `workDir` and `outdir` parameters in `configs/nextflow.base.config`
* Edit the `runOptions` in `configs/nextflow.base.config` in order to mount directories to your run

When analysing Nanopore data:
* Edit the `ext.args` for Flye: specify genome size for the organism of interest with flag `--genome-size`
* Edit the `ext.seqmethod`for Flye depending on the input data
* Edit the `ext.args` for Medaka: specify the model with flag `-m`. Currently it is set to `r941_min_sup_g507`, but one should always set it based on how the data was produced. More about choosing the right model can be found [here](https://github.com/nanoporetech/medaka#models).

### Test data
Source: `assets/test_data/samplelist.csv`

* Edit the read1 and read2 columns in `assets/test_data/samplelist.csv`

## Setting up temp directories

Source: `~/.bashrc`

* Add the export line to `~/.bashrc`
* Change `SINGULARITY_TMPDIR` to `APPTAINER_TMPDIR` if you are using apptainer

```bash
export SINGULARITY_TMPDIR="/tmp" #or equivalent filepath to tmp dir
```

## Fetching databases

### Choose database

Choose between Kraken DB (64GB [Highly recommended]) or MiniKraken DB (8GB).  Or customize [your own](https://benlangmead.github.io/aws-indexes/k2).

### Download Kraken database

```bash
wget -O /path/to/kraken_db/krakenstd.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenstd.tar.gz
```

### Download MiniKraken database

```bash
wget -O /path/to/kraken_db/krakenmini.tar.gz https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230314.tar.gz
tar -xf /path/to/kraken_db/krakenmini.tar.gz
```

## Updating databases

### Update MLST database

```bash
bash /path/to/jasen/assets/mlst_db/update_mlst_db.sh
```

## Create personalised TBProfiler database

### Install jasentool

```bash
git clone git@github.com:ryanjameskennedy/jasentool.git && cd jasentool
pip install .
```

### Bgzip and index gms TBProfiler db

```bash
bgzip -c converged_who_fohm_tbdb.bed > /path/to/jasen/assets/tbprofiler_dbs/bed/converged_who_fohm_tbdb.bed.gz
tabix -p bed /path/to/jasen/assets/tbprofiler_dbs/bed/converged_who_fohm_tbdb.bed.gz
```
